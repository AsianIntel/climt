module SlabSurface {
    use Components;
    use Map;
    use Properties;
    use marker;
    use State;
    use DataArray;

    type FluxMarker = UnitMarker(0, 1, -3, 0, 0, 0, 0);
    type dimensionless = UnitMarker(0, 0, 0, 0, 0, 0, 0);

    class SlabSurface: TendencyComponent {
        proc init() {
            var input_properties = new map(string, shared AbstractProperties);

            input_properties.add("downwelling_longwave_flux_in_air", new shared Properties({"*", "interface_levels"}, new FluxMarker(1, 0, "W m^-2")));
            input_properties.add("downwelling_shortwave_flux_in_air", new shared Properties({"*", "interface_levels"}, new FluxMarker(1, 0, "W m^-2")));
            input_properties.add("upwelling_longwave_flux_in_air", new shared Properties({"*", "interface_levels"}, new FluxMarker(1, 0, "W m^-2")));
            input_properties.add("upwelling_shortwave_flux_in_air", new shared Properties({"*", "interface_levels"}, new FluxMarker(1, 0, "W m^-2")));
            input_properties.add("surface_upward_latent_heat_flux", new shared Properties({"*"}, new FluxMarker(1, 0, "W m^-2")));
            input_properties.add("surface_temperature", new shared Properties({"*"}, new UnitMarker(0, 0, 0, 0, 1, 0, 0, 1, 0, "K")));
            input_properties.add("surface_upward_sensible_heat_flux", new shared Properties({"*"}, new FluxMarker(1, 0, "W m^-2")));
            input_properties.add("surface_thermal_capacity", new shared Properties({"*"}, new UnitMarker(2, 0, -2, 0, -1, 0, 0, 1, 0, "J kg^-1 degK^-1")));
            input_properties.add("surface_material_density", new shared Properties({"*"}, new UnitMarker(-3, 1, 0, 0, 0, 0, 0, 1, 0, "kg m^-3")));
            input_properties.add("upward_heat_flux_at_ground_level_in_soil", new shared Properties({"*"}, new FluxMarker(1, 0, "W m^-2")));
            input_properties.add("heat_flux_into_sea_water_due_to_sea_ice", new shared Properties({"*"}, new FluxMarker(1, 0, "W m^-2")));
            input_properties.add("area_type", new shared Properties({"*"}, new dimensionless(1, 0, "")));
            input_properties.add("soil_layer_thickness", new shared Properties({"*"}, new UnitMarker(1, 0, 0, 0, 0, 0, 0, 1, 0, "m")));
            input_properties.add("ocean_mixed_layer_thickness", new shared Properties({"*"}, new UnitMarker(1, 0, 0, 0, 0, 0, 0, 1, 0, "m")));
            input_properties.add("heat_capacity_of_soil", new shared Properties({"*"}, new UnitMarker(2, 0, -2, 0, -1, 0, 0, 1, 0, "J kg^-1 degK^-1")));
            input_properties.add("sea_water_density", new shared Properties({"*"}, new UnitMarker(-3, 1, 0, 0, 0, 0, 0, 1, 0, "kg m^-3")));

            var tendency_properties = new map(string, shared AbstractProperties);
            tendency_properties.add("surface_temperature", new shared Properties({"*"}, new UnitMarker(0, 0, -1, 0, 1, 0, 0, 1, 0, "K s^-1")));

            var diagnostic_properties = new map(string, shared AbstractProperties);
            diagnostic_properties.add("depth_of_slab_surface", new shared Properties({"*"}, new UnitMarker(1, 0, 0, 0, 0, 0, 0, 1, 0, "m")));

            super.init(input_properties, tendency_properties, diagnostic_properties);
        }

        proc array_call(state: State) {
            var diagnostics = initialize_arrays_with_properties(this.diagnostic_properties, state, this.input_properties);
            var tendencies = initialize_arrays_with_properties(this.tendency_properties, state, this.input_properties);

            var downwelling_shortwave_flux_in_air = try! state.getValue("downwelling_shortwave_flux_in_air"): DataArray2(real, FluxMarker);
            var downwelling_longwave_flux_in_air = try! state.getValue("downwelling_longwave_flux_in_air"): DataArray2(real, FluxMarker);
            var upwelling_shortwave_flux_in_air = try! state.getValue("upwelling_shortwave_flux_in_air"): DataArray2(real, FluxMarker);
            var upwelling_longwave_flux_in_air = try! state.getValue("upwelling_longwave_flux_in_air"): DataArray2(real, FluxMarker);
            var surface_upward_sensible_heat_flux = try! state.getValue("surface_upward_sensible_heat_flux"): DataArray1(real, FluxMarker);
            var surface_upward_latent_heat_flux = try! state.getValue("surface_upward_latent_heat_flux"): DataArray1(real, FluxMarker);

            var net_heat_flux = (
                downwelling_shortwave_flux_in_air.arr(.., 0) +
                downwelling_longwave_flux_in_air.arr(.., 0) -
                upwelling_shortwave_flux_in_air.arr(.., 0) -
                upwelling_longwave_flux_in_air.arr(.., 0) -
                surface_upward_sensible_heat_flux.arr -
                surface_upward_latent_heat_flux.arr
            );
            
            var area_type = try! state.getValue("area_type"): DataArray1(int, dimensionless);
            var land_mask = area_type.arr[0] == 0 /* land */ || area_type.arr[0] == 1 /* land_ice */: int;
            var sea_mask = area_type.arr[0] == 2 /* sea */ || area_type.arr[0] == 3 /* sea_ice */: int;
            var land_ice_mask = area_type.arr[0] == 1 /* land_ice */: int;
            var sea_ice_mask = area_type.arr[0] == 3 /* sea_ice */: int;

            var upward_heat_flux_at_ground_level_in_soil = try! state.getValue("upward_heat_flux_at_ground_level_in_soil"): DataArray1(real, FluxMarker);
            net_heat_flux[land_ice_mask] = -upward_heat_flux_at_ground_level_in_soil.arr[0];
            var heat_flux_into_sea_water_due_to_sea_ice = try! state.getValue("heat_flux_into_sea_water_due_to_sea_ice"): DataArray1(real, FluxMarker);
            net_heat_flux[sea_ice_mask] = heat_flux_into_sea_water_due_to_sea_ice.arr[sea_ice_mask];
            
            var surface_material_density = try! state.getValue("surface_material_density"): DataArray1(real, UnitMarker(-3, 1, 0, 0, 0, 0, 0));
            var sea_water_density = try! state.getValue("sea_water_density"): DataArray1(real, UnitMarker(-3, 1, 0, 0, 0, 0, 0));
            surface_material_density.arr[sea_mask] = sea_water_density.arr[sea_mask];

            var surface_thermal_capacity = try! state.getValue("surface_thermal_capacity"): DataArray1(real, UnitMarker(2, 0, -2, 0, -1, 0, 0));
            var heat_capacity_of_soil = try! state.getValue("heat_capacity_of_soil"): DataArray1(real, UnitMarker(2, 0, -2, 0, -1, 0, 0));
            surface_thermal_capacity.arr[land_mask] = heat_capacity_of_soil.arr[land_mask];

            var depth_of_slab_surface = try! diagnostics.getValue("depth_of_slab_surface"): DataArray1(real, UnitMarker(1, 0, 0, 0, 0, 0, 0));
            var ocean_mixed_layer_thickness = try! state.getValue("ocean_mixed_layer_thickness"): DataArray1(real, UnitMarker(1, 0, 0, 0, 0, 0, 0));
            var soil_layer_thickness = try! state.getValue("soil_layer_thickness"): DataArray1(real, UnitMarker(1, 0, 0, 0, 0, 0, 0));
            depth_of_slab_surface.arr[sea_mask] = ocean_mixed_layer_thickness.arr[sea_mask];
            depth_of_slab_surface.arr[land_mask] = soil_layer_thickness.arr[land_mask];

            var mass_surface_slab = surface_material_density.arr * depth_of_slab_surface.arr;
            var heat_capacity_surface = mass_surface_slab * surface_thermal_capacity.arr;

            var surface_temperature = try! tendencies.getValue("surface_temperature"): DataArray1(real, UnitMarker(0, 0, -1, 0, 1, 0, 0));
            surface_temperature.arr(..) = net_heat_flux / heat_capacity_surface;
            surface_temperature.arr(land_ice_mask) = 0.0;
            surface_temperature.arr(sea_ice_mask) = 0.0;

            return (tendencies, diagnostics);
        }
    }
}