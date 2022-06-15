module HybridSigmaPressure {
    use Components;
    use Map;
    use Properties;
    use State;
    use marker;
    use DataArray;
    use LinearAlgebra;

    class HybridSigmaPressureDiagnosticComponent: DiagnosticComponent {
        proc init() {
            var input_properties = new map(string, shared Properties);
            input_properties.add("atmosphere_hybrid_sigma_pressure_a_coordinate_on_interface_levels", new shared Properties({"interface_levels", "*"}, new UnitMarker(0, 0, 0, 0, 0, 0, 0, 1, 0, "")));
            input_properties.add("atmosphere_hybrid_sigma_pressure_b_coordinate_on_interface_levels", new shared Properties({"interface_levels", "*"}, new UnitMarker(0, 0, 0, 0, 0, 0, 0, 1, 0, "")));
            input_properties.add("surface_air_pressure", new shared Properties({"*"}, new UnitMarker(-1, 1, -2, 0, 1, 0, 0, 1, 0, "Pa")));

            var diagnostic_properties = new map(string, shared Properties);
            diagnostic_properties.add("air_pressure", new shared Properties({"mid_levels", "*"}, new UnitMarker(-1, 1, -2, 0, 1, 0, 0, 1, 0, "Pa")));
            diagnostic_properties.add("air_pressure_on_interface_levels", new shared Properties({"interface_levels", "*"}, new UnitMarker(-1, 1, -2, 0, 1, 0, 0, 1, 0, "Pa")));
            super.init(input_properties, diagnostic_properties);
        }

        proc array_call(state: State) {
            var model_top_pressure = 20.0;

            var a_coord = try! state.getValue("atmosphere_hybrid_sigma_pressure_a_coordinate_on_interface_levels"): DataArray2(real);
            var b_coord = try! state.getValue("atmosphere_hybrid_sigma_pressure_b_coordinate_on_interface_levels"): DataArray2(real);
            var surface_air_pressure = try! state.getValue("surface_air_pressure"): DataArray2(real);

            var sap = surface_air_pressure.arr - model_top_pressure;
            var b = b_coord.arr.dot(sap);
            var p_interface = a_coord.arr + b_coord.arr.dot(sap);

            var p_dims = p_interface.shape;
            var delta_p: [0..p_dims[0]-2, 0..p_dims[1]-1] real;
            for i in 0..p_dims[0]-2 {
                delta_p(i, ..) = p_interface(i+1, ..) - p_interface(i, ..);
            }

            var gas_constant_of_dry_air = 287.0;
            var heat_capacity_of_dry_air_at_constant_pressure = 1004.64;
            var rk = gas_constant_of_dry_air / heat_capacity_of_dry_air_at_constant_pressure;

            var p: [0..p_dims[0]-2, 0..p_dims[1]-1] real;
            for i in 0..p_dims[0]-2 {
                p(i, ..) = (
                    (p_interface(i+1, ..) ** (rk + 1) - p_interface(i, ..) ** (rk + 1)) / ((rk + 1) * delta_p(i, ..))
                ) ** (1/rk);
            }
            // writeln(p);

            state.add("air_pressure", new shared DataArray2(p, a_coord.dimensions));
            state.add("air_pressure_on_interface_levels", new shared DataArray2(p_interface, a_coord.dimensions));
        }
    }
}