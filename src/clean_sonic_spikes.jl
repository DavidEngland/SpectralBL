"""
    clean_sonic_spikes(signal::Vector{Float64}; threshold::Float64=3.5)

Applies a robust Median Absolute Deviation (MAD) pass to strip spikes out of
raw boundary-layer time series data before calculating 10-minute statistics.
"""
function clean_sonic_spikes(signal::Vector{Float64}; threshold::Float64=3.5)
    med = median(signal)
    # Calculate MAD: median(|x_i - median(x)|)
    abs_dev = map(x -> abs(x - med), signal)
    mad_val = median(abs_dev)

    # Prevent divide-by-zero on perfectly flat signal lines
    if mad_val == 0.0
        return signal
    end

    # 0.6745 normalizes MAD to make it comparable to standard deviation
    modified_z_scores = map(x -> (0.6745 * (x - med)) / mad_val, signal)

    # Generate clean output vector, substituting extreme spikes with NaN
    # to maintain temporal spacing for subsequent spectral operations
    cleaned_signal = copy(signal)
    spike_indices = findall(z -> abs(z) > threshold, modified_z_scores)

    if !isempty(spike_indices)
        cleaned_signal[spike_indices] .= NaN
        # Optional tracing for data auditing logs
        # println("  [FILTER] Suppressed $(length(spike_indices)) extreme signal spikes.")
    end

    return cleaned_signal
end