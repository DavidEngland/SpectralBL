using LinearAlgebra
using Statistics
using DataFrames
using CSV
using Plots

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using UnifiedManifold

function chebyshev_basis_matrix(xi::Vector{Float64}, N::Int)
    B = zeros(Float64, length(xi), N + 1)
    for i in eachindex(xi)
        θ = acos(clamp(xi[i], -1.0, 1.0))
        for n in 0:N
            B[i, n + 1] = cos(n * θ)
        end
    end
    return B
end

function spectral_damp!(u_state::Vector{Float64}, v_state::Vector{Float64}, B::Matrix{Float64}, psi_T::Vector{Float64}, damping_rate::Float64)
    a_u = B \ u_state
    a_v = B \ v_state

    for i in eachindex(a_u)
        # Damping envelope targets only high-frequency turbulence-window modes.
        damp = exp(-damping_rate * psi_T[i])
        a_u[i] *= damp
        a_v[i] *= damp
    end

    u_state .= B * a_u
    v_state .= B * a_v
end

function rk4_step!(u_state::Vector{Float64}, v_state::Vector{Float64}, D2::Matrix{Float64}, c_wave::Float64, dt::Float64)
    k1_u = v_state
    k1_v = c_wave^2 .* (D2 * u_state)

    k2_u = v_state .+ 0.5 * dt .* k1_v
    k2_v = c_wave^2 .* (D2 * (u_state .+ 0.5 * dt .* k1_u))

    k3_u = v_state .+ 0.5 * dt .* k2_v
    k3_v = c_wave^2 .* (D2 * (u_state .+ 0.5 * dt .* k2_u))

    k4_u = v_state .+ dt .* k3_v
    k4_v = c_wave^2 .* (D2 * (u_state .+ dt .* k3_u))

    u_state .+= (dt / 6.0) .* (k1_u .+ 2.0 .* k2_u .+ 2.0 .* k3_u .+ k4_u)
    v_state .+= (dt / 6.0) .* (k1_v .+ 2.0 .* k2_v .+ 2.0 .* k3_v .+ k4_v)
end

function stable_timestep(D2::Matrix{Float64}, c_wave::Float64, dt_nominal::Float64)
    λ = eigvals(c_wave^2 .* D2)
    ωmax = sqrt(maximum(abs.(λ)) + 1e-12)
    return min(dt_nominal, 0.2 / ωmax)
end

function run_scenario!(u_initial::Vector{Float64}, v_initial::Vector{Float64}, D2::Matrix{Float64}, B::Matrix{Float64}, ws, c_wave::Float64, dt::Float64, n_steps::Int; apply_sponge::Bool, damping_rate::Float64, sample_stride::Int)
    u_state = copy(u_initial)
    v_state = copy(v_initial)
    top_energy = Float64[]

    for step in 1:n_steps
        rk4_step!(u_state, v_state, D2, c_wave, dt)

        if !all(isfinite, u_state) || !all(isfinite, v_state)
            error("Non-finite state encountered at step $step. Reduce dt or damping strength.")
        end

        if apply_sponge
            spectral_damp!(u_state, v_state, B, ws.psi_T, damping_rate)
        end

        if step % sample_stride == 0
            # ws.z_atm is ordered from top to bottom in this workspace.
            push!(top_energy, 0.5 * (u_state[1]^2 + v_state[1]^2))
        end
    end

    return top_energy
end

function test_gravity_wave_sponge_layer()
    ws = UnifiedManifoldWorkspace(32, 1.5, 50.0, 0.05, K_q = 72)
    N = ws.N
    B = chebyshev_basis_matrix(Vector{Float64}(ws.xi_target), N)

    c_wave = 2.0
    z_center = 35.0
    sigma_z = 2.0
    wavelength = 8.0

    u_initial = [exp(-0.5 * ((z - z_center) / sigma_z)^2) * sin(2.0 * pi * (z - z_center) / wavelength) for z in ws.z_atm]
    u_initial = Vector{Float64}(u_initial)

    v_initial = -c_wave .* (Matrix{Float64}(ws.Dz_atm) * u_initial)

    D2_raw = Matrix{Float64}(ws.Dz_atm) * Matrix{Float64}(ws.Dz_atm)
    D2_sym = 0.5 .* (D2_raw .+ transpose(D2_raw))
    eig = eigen(Symmetric(D2_sym))
    λ_stable = min.(eig.values, 0.0)
    D2 = Matrix(eig.vectors * Diagonal(λ_stable) * transpose(eig.vectors))

    dt_nominal = 0.02
    dt = stable_timestep(D2, c_wave, dt_nominal)
    total_time = 100.0
    n_steps = max(2000, Int(round(total_time / dt)))
    sample_stride = max(1, n_steps ÷ 100)
    damping_rate = 0.05

    println("Using stable dt = ", round(dt, sigdigits = 4), " with ", n_steps, " steps")

    energy_with = run_scenario!(u_initial, v_initial, D2, B, ws, c_wave, dt, n_steps;
                                apply_sponge = true, damping_rate = damping_rate, sample_stride = sample_stride)
    energy_without = run_scenario!(u_initial, v_initial, D2, B, ws, c_wave, dt, n_steps;
                                   apply_sponge = false, damping_rate = 0.0, sample_stride = sample_stride)

    n_tail = min(20, length(energy_with), length(energy_without))
    E_masked = mean(energy_with[end - n_tail + 1:end])
    E_unmasked = mean(energy_without[end - n_tail + 1:end])
    suppression_ratio = E_unmasked / (E_masked + 1e-12)

    println("=== SPONGE LAYER EXPERIMENT RESULTS ===")
    println("Mean Boundary Energy (WITHOUT Sponge): ", round(E_unmasked, digits = 6))
    println("Mean Boundary Energy (WITH Sponge):    ", round(E_masked, digits = 6))
    println("Reflection Suppression Performance:    ", round(suppression_ratio, digits = 2), "x Ratio")

    time_axis = collect(1:length(energy_with))
    p = plot(time_axis, energy_with, label = "WITH Sponge (Modal Masking)", linewidth = 2,
             xlabel = "Sample Window (Time)", ylabel = "Energy at Upper Boundary",
             title = "Sponge Layer Validation Experiment", legend = :topright,
             left_margin = 14Plots.mm, bottom_margin = 8Plots.mm, size = (1100, 550))
    plot!(p, time_axis, energy_without, label = "WITHOUT Sponge (Free Reflection)", linewidth = 2)

    mkpath("data")
    savefig(p, "data/wave_reflection_test.png")
    savefig(p, "data/wave_reflection_test.pdf")

    diag = DataFrame(
        sample = time_axis,
        energy_top_with_sponge = energy_with,
        energy_top_without_sponge = energy_without,
    )
    CSV.write("data/wave_reflection_metrics.csv", diag)

    println("✓ Plot saved to: data/wave_reflection_test.png")
    println("✓ Plot saved to: data/wave_reflection_test.pdf")
    println("✓ Metrics saved to: data/wave_reflection_metrics.csv")

    return suppression_ratio
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_gravity_wave_sponge_layer()
end
