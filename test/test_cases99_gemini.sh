julia --project -e '
    include("src/Cases99.gemini.jl");
    using .UnifiedManifold;
    ws = UnifiedManifoldWorkspace(32, 1.5, 429.5, 2.5);
    println("Active Map Type: ", typeof(ws.map));
    @assert typeof(ws.map) <: CoordinateMap "Validation Failed: Custom mapping type out of bounds";
    println("✅ [SUCCESS] Workspace unifies perfectly under TanhMap default.");
'