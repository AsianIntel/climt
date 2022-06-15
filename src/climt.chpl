/* Documentation for climt */
module climt {
  use Map;
  import Math;
  use RRTMGLongwave;
  use Properties;
  use Components;
  import utilities;
  use State;
  use DataArray;
  use marker;
  use HybridSigmaPressure;
  use SlabSurface;

  writeln("New library: climt");

  proc get_atmospheric_grid(name: string, value: real, state: State, interface_: bool = false, horizontal: bool = false) {
    var latitude = try! state.getValue("latitude"): DataArray2(real);
    var shape = latitude.arr.shape;

    var hybrid_sigma_pressure = try! state.getValue("atmosphere_hybrid_sigma_pressure_a_coordinate_on_interface_levels"): DataArray2(real);

    if horizontal {
      var arr: [0..shape[1]-1, 0..shape[0]-1] real = value;
      state.add(name, new shared DataArray2(arr, latitude.dimensions));
    } else if interface_ {
      var arr: [0..hybrid_sigma_pressure.arr.shape[0]-1, 0..shape[1]-1, 0..shape[0]-1] real = value;
      var hsd = latitude.dimensions.these();
      var dom = {"interface_levels", hsd[0], hsd[1]};
      state.add(name, new shared DataArray3(arr, dom));
    } else {
      var arr: [0..hybrid_sigma_pressure.arr.shape[0]-2, 0..shape[1]-1, 0..shape[0]-1] real = value;
      var hsd = latitude.dimensions.these();
      var dom = {"mid_levels", hsd[0], hsd[1]};
      state.add(name, new shared DataArray3(arr, dom));
    }
  }

  proc get_surface_grid(name: string, value, state: State) {
    var latitude = try! state.getValue("latitude"): DataArray2(real);
    var shape = latitude.arr.shape;

    var arr: [0..shape[1]-1, 0..shape[0]-1] value.type = value;
    state.add(name, new shared DataArray2(arr, latitude.dimensions));
  }

  proc get_land_grid(name: string, value: real, state: State, horizontal: bool = false) {
    var latitude = try! state.getValue("latitude"): DataArray2(real);
    var shape = latitude.arr.shape;

    if horizontal {
      var arr: [0..shape[1]-1, 0..shape[0]-1] value.type = value;
      state.add(name, new shared DataArray2(arr, latitude.dimensions));
    }
  }

  proc get_ocean_grid(name: string, value: real, state: State, horizontal: bool = false) {
    var latitude = try! state.getValue("latitude"): DataArray2(real);
    var shape = latitude.arr.shape;

    if horizontal {
      var arr: [0..shape[1]-1, 0..shape[0]-1] real = value;
      state.add(name, new shared DataArray2(arr, latitude.dimensions));
    }
  }

  proc get_ice_grid(name: string, value: real, state: State, interface_: bool = false, horizontal: bool = false) {
    var latitude = try! state.getValue("latitude"): DataArray2(real);
    var shape = latitude.arr.shape;

    var height_on_ice = try! state.getValue("height_on_ice_interface_levels"): DataArray1(real);

    if horizontal {
      var arr: [0..shape[1]-1, 0..shape[0]-1] real = value;
      state.add(name, new shared DataArray2(arr, latitude.dimensions));
    } else if interface_ {
      var arr: [0..height_on_ice.arr.shape[0]-1, 0..shape[1]-1, 0..shape[0]-1] real = value;
      var hsd = latitude.dimensions.these();
      var dom = {"ice_interface_levels", hsd[0], hsd[1]};
      state.add(name, new shared DataArray3(arr, dom));
    } else {
      var arr: [0..height_on_ice.arr.shape[0]-2, 0..shape[1]-1, 0..shape[0]-1] real = value;
      var hsd = latitude.dimensions.these();
      var dom = {"ice_mid_levels", hsd[0], hsd[1]};
      state.add(name, new shared DataArray3(arr, dom));
    }
  }

  proc default_values(name: string, state: State) {
    select name {
      when "air_temperature" {
        get_atmospheric_grid(name, 290.0, state, false, false);
      } 
      when "northward_wind" {
        get_atmospheric_grid(name, 0.0, state, false, false);
      } 
      when "eastward_wind" {
        get_atmospheric_grid(name, 0.0, state, false, false);
      } 
      when "divergence_of_wind" {
        get_atmospheric_grid(name, 0.0, state, false, false);
      } 
      when "atmospheric_relative_vorticity" {
        get_atmospheric_grid(name, 0.0, state, false, false);
      } 
      when "specific_humidity" {
        get_atmospheric_grid(name, 0.0, state, false, false);
      } 
      when "mole_fraction_of_carbon_dioxide_in_air" {
        get_atmospheric_grid(name, 330e-6, state, false, false);
      } 
      when "mole_fraction_of_methane_in_air" {
        get_atmospheric_grid(name, 0.0, state, false, false);
      } 
      when "mole_fraction_of_nitrous_oxide_in_air" {
        get_atmospheric_grid(name, 0.0, state, false, false);
      } 
      when "mole_fraction_of_oxygen_in_air" {
        get_atmospheric_grid(name, 0.21, state, false, false);
      } 
      when "mole_fraction_of_nitrogen_in_air" {
        get_atmospheric_grid(name, 0.78, state, false, false);
      } 
      when "mole_fraction_of_hydrogen_in_air" {
        get_atmospheric_grid(name, 500e-9, state, false, false);
      } 
      when "mole_fraction_of_cfc11_in_air" {
        get_atmospheric_grid(name, 0.0, state, false, false);
      } 
      when "mole_fraction_of_cfc12_in_air" {
        get_atmospheric_grid(name, 0.0, state, false, false);
      } 
      when "mole_fraction_of_cfc22_in_air" {
        get_atmospheric_grid(name, 0.0, state, false, false);
      } 
      when "mole_fraction_of_carbon_tetrachloride_in_air" {
        get_atmospheric_grid(name, 0.0, state, false, false);
      } 
      when "cloud_area_fraction_in_atmosphere_layer" {
        get_atmospheric_grid(name, 0.0, state, false, false);
      } 
      when "mass_content_of_cloud_ice_in_atmosphere_layer" {
        get_atmospheric_grid(name, 0.0, state, false, false);
      } 
      when "mass_content_of_cloud_liquid_water_in_atmosphere_layer" {
        get_atmospheric_grid(name, 0.0, state, false, false);
      } 
      when "cloud_ice_particle_size" {
        get_atmospheric_grid(name, 20.0, state, false, false);
      } 
      when "cloud_water_droplet_radius" {
        get_atmospheric_grid(name, 10.0, state, false, false);
      } 
      when "cloud_base_mass_flux" {
        get_atmospheric_grid(name, 0.0, state, false, true);
      } 
      when "zenith_angle" {
        get_atmospheric_grid(name, 0.0, state, false, true);
      } 
      when "downwelling_shortwave_flux_in_air" {
        get_atmospheric_grid(name, 0.0, state, true, false);
      } 
      when "downwelling_longwave_flux_in_air" {
        get_atmospheric_grid(name, 0.0, state, true, false);
      } 
      when "upwelling_shortwave_flux_in_air" {
        get_atmospheric_grid(name, 0.0, state, true, false);
      } 
      when "upwelling_longwave_flux_in_air" {
        get_atmospheric_grid(name, 0.0, state, true, false);
      } 
      when "surface_specific_humidity" {
        get_surface_grid(name, 0.0, state);
      } 
      when "surface_temperature" {
        get_surface_grid(name, 300.0, state);
      } 
      when "soil_surface_temperature" {
        get_surface_grid(name, 300.0, state);
      } 
      when "surface_geopotential" {
        get_surface_grid(name, 0.0, state);
      } 
      when "surface_thermal_capacity" {
        get_surface_grid(name, 4.1813e3, state);
      } 
      when "depth_of_slab_surface" {
        get_surface_grid(name, 50.0, state);
      } 
      when "surface_material_density" {
        get_surface_grid(name, 1000.0, state);
      } 
      when "surface_albedo_for_direct_shortwave" {
        get_surface_grid(name, 0.06, state);
      } 
      when "surface_albedo_for_diffuse_shortwave" {
        get_surface_grid(name, 0.06, state);
      } 
      when "surface_albedo_for_direct_near_infrared" {
        get_surface_grid(name, 0.06, state);
      } 
      when "surface_albedo_for_diffuse_near_infrared" {
        get_surface_grid(name, 0.06, state);
      } 
      when "surface_roughness_length" {
        get_surface_grid(name, 0.0002, state);
      } 
      when "surface_drag_coefficient_for_heat_in_air" {
        get_surface_grid(name, 0.0012, state);
      } 
      when "surface_drag_coefficient_for_momentum_in_air" {
        get_surface_grid(name, 0.0012, state);
      } 
      when "area_type" {
        get_surface_grid(name, 0: int, state);
      } 
      when "surface_upward_sensible_heat_flux" {
        get_surface_grid(name, 0.0, state);
      } 
      when"surface_upward_latent_heat_flux" {
        get_surface_grid(name, 0.0, state);
      } 
      when "soil_type" {
        get_land_grid(name, 0: int, state, true);
      } 
      when "soil_temperature" {
        get_land_grid(name, 274.0, state, true);
      } 
      when "soil_layer_thickness" {
        get_land_grid(name, 50.0, state, true);
      } 
      when "upward_heat_flux_at_ground_level_in_soil" {
        get_land_grid(name, 0.0, state, true);
      } 
      when "heat_capacity_of_soil" {
        get_land_grid(name, 2000.0, state, true);
      } 
      when "sea_water_density" {
        get_ocean_grid(name, 1.029e3, state, true);
      } 
      when "sea_surface_temperature" {
        get_ocean_grid(name, 300.0, state, true);
      } 
      when "ocean_mixed_layer_thickness" {
        get_ocean_grid(name, 50.0, state, true);
      } 
      when "snow_and_ice_temperature" {
        get_ice_grid(name, 270.0, state, true, false);
      } 
      when "heat_flux_into_sea_water_due_to_sea_ice" {
        get_ice_grid(name, 0.0, state, false, true);
      }
      when "land_ice_thickness" {
        get_ice_grid(name, 0.0, state, false, true);
      } 
      when "sea_ice_thickness" {
        get_ice_grid(name, 0.0, state, false, true);
      } 
      when "surface_snow_thickness" {
        get_ice_grid(name, 0.0, state, false, true);
      } 
      when "solar_cycle_fraction" {
        var arr = [0.0];
        var dom: domain(string);
        state.add(name, new shared DataArray1(arr, dom));
      } 
      when "flux_adjustment_for_earth_sun_distance" {
        var arr = [1.0];
        var dom: domain(string);
        state.add(name, new shared DataArray1(arr, dom));
      } 
      when "lwe_thickness_of_soil_moisture_content" {
        get_surface_grid(name, 0.0, state);
      } 
      when "convective_precipitation_rate" {
        get_surface_grid(name, 0.0, state);
      } 
      when "stratiform_precipitation_rate" {
        get_surface_grid(name, 0.0, state);
      }
    }
  }

  proc get_grid(
    nx: int = 1, ny: int = 1, nz: int = 28, n_ice_interface_levels: int = 10, p_surf_in_Pa: real = 1.0132e5, 
    p_toa_in_Pa: real = 20.0, proportion_sigma_levels: real = 0.1, proportion_isobaric_level: real = 0.25,
    x_name: string = "lon", y_name: string = "lat", latitude_grid: string = "regular"
  ) {
    var return_state = get_hybrid_sigma_pressure_levels(nz + 1, p_surf_in_Pa, p_toa_in_Pa, proportion_isobaric_level, proportion_sigma_levels);

    var surface_arr: [0..ny-1, 0..nx-1] real = 1.0 * p_surf_in_Pa;
    return_state.add("surface_air_pressure", new shared DataArray2(surface_arr, {y_name, x_name}));

    var component = new HybridSigmaPressureDiagnosticComponent();
    component.array_call(return_state);

    var two_dim_lons: [0..ny-1, 0..nx-1] real = 1.0;
    var space = utilities.linspace(0.0, 360.0, nx*2, false);
    for i in 0..ny-1 {
      var j = 0;
      while (j <= nx-1) {
        two_dim_lons(i, j) = space[j*2]; 
        j += 1;
      }
    }
    return_state.add("longitude", new shared DataArray2(two_dim_lons, {y_name, x_name}));

    var two_dim_lats: [0..ny-1, 0..nx-1] real = 1.0;
    if latitude_grid == "regular" {
      var space_lats = utilities.linspace(-90.0, 90.0, ny * 2 + 1, true);
      for i in 0..ny-1 {
        var j = 0;
        while (j <= nx - 1) {
          two_dim_lats(i, j) = space_lats[j * 2];
          j += 1;
        }
      } 
      return_state.add("latitude", new shared DataArray2(two_dim_lats, {y_name, x_name}));
    } else if latitude_grid == "gaussian" {
      // TODO
    }

    var ice_interface_levels: [0..n_ice_interface_levels-1] real = 0.0;
    return_state.add("height_on_ice_interface_levels", new shared DataArray1(ice_interface_levels, {"ice_interface_levels"}));

    return return_state;
  }

  proc get_hybrid_sigma_pressure_levels(
    num_levels: int = 28, reference_pressure: real = 1e5, model_top_pressure: real = 20.0,
    proportion_isobaric_level: real = 0.25, proportion_sigma_levels: real = 0.1
  ) {
    var thickness_dist = utilities.sin(utilities.linspace(0.1, Math.pi-0.1, num_levels-1));
    thickness_dist /= utilities.sum(thickness_dist);
    thickness_dist *= (reference_pressure - model_top_pressure);

    var pressure_levels: [0..num_levels-1] real = 0.0;
    pressure_levels[0] = model_top_pressure;

    var sum = thickness_dist[0];
    for i in 1..num_levels-2 {
      pressure_levels[i] = model_top_pressure + sum;
      sum += thickness_dist[i];
    }
    pressure_levels[num_levels-1] = model_top_pressure + sum;

    var sigma_interface = (pressure_levels - model_top_pressure)/(reference_pressure - model_top_pressure);

    var ak: [0..num_levels-1] real = 0.0;
    var bk: [0..num_levels-1] real = 0.0;

    var num_isobaric_levels = (proportion_isobaric_level * num_levels): int;
    var num_sigma_levels = (proportion_sigma_levels * num_levels): int;

    for i in 0..num_isobaric_levels {
      ak[i] = pressure_levels[i];
    }

    var isobaric_sigma_level = sigma_interface[num_isobaric_levels - 1];

    for level in num_isobaric_levels..(num_levels-num_sigma_levels-1) {
      var sigma_value = sigma_interface[level];
      var b_level = (sigma_value - isobaric_sigma_level) / (1 - isobaric_sigma_level);
      var r_level = get_exponent_for_sigma(b_level, num_sigma_levels);
      bk[level] = b_level ** r_level;
      ak[level] = model_top_pressure + (sigma_value - bk[level]) * (reference_pressure - model_top_pressure);
    }

    for level in (num_levels - num_sigma_levels)..num_levels-1 {
      var sigma_value = sigma_interface[level];
      bk[level] = (sigma_interface[level] - isobaric_sigma_level) / (1 - isobaric_sigma_level);
      ak[level] = model_top_pressure + (sigma_value - bk[level]) * (reference_pressure - model_top_pressure);
    }

    ak.reverse();
    bk.reverse();

    var ak2 = reshape(ak, {0..num_levels-1, 0..0});
    var bk2 = reshape(bk, {0..num_levels-1, 0..0});

    var state = new State();
    state.add("atmosphere_hybrid_sigma_pressure_a_coordinate_on_interface_levels", new shared DataArray2(ak2, {"interface_levels", "*"}));
    state.add("atmosphere_hybrid_sigma_pressure_b_coordinate_on_interface_levels", new shared DataArray2(bk2, {"interface_levels", "*"}));

    return state;
  }

  proc get_exponent_for_sigma(b_half, num_sigma_levels) {
    var r_p = 2.2;
    var r_sigma = 1.0;
    var S = 5;

    if num_sigma_levels > 0 {
      r_sigma = 1.0;
    } else {
      r_sigma = 1.35;
    }

    return r_p + (r_sigma - r_p) * Math.atan(S * b_half) / Math.atan(S);
  }

  proc get_default_state(in components, n_ice_interface_levels: int = 30) {
    var grid_state = get_grid(n_ice_interface_levels=n_ice_interface_levels);

    var properties = new map(string, shared Properties);
    for component in components {
      for property in component.input_properties.items() {
        properties.add(property[0], property[1]);
      }
    }
    
    for name in properties.keys() {
      if !grid_state.data_map.contains(name) {
        default_values(name, grid_state);
      }
    }

    return grid_state;
  }

  proc get_dim_lengths_from_raw_input(raw_input: State, input_properties: map(string, shared Properties)) throws {
    var dim_lengths = new map(string, int);

    for (name, property) in input_properties.items() {
        var i = 0;
        var ada = raw_input.getValue(name);
        for dim_name in ada.getDims() {
          dim_lengths.add(dim_name, ada.getDomShape()[i]);
          i += 1;
        }   
    }

    return dim_lengths;
  }

  proc initialize_arrays_with_properties(output_properties: map(string, shared Properties), input_state: State, input_properties: map(string, shared Properties)) {
    var dim_lengths = get_dim_lengths_from_raw_input(input_state, input_properties);

    var out_dict = new map(string, shared AbstractDataArray);
    for (name, property) in output_properties.items() {
      var out_shape: [property.dims.size] int;
      for i in 0..property.dims.size-1 {
        var dim = property.dims[i];
        out_shape[i] = dim_lengths[dim];
        if (out_shape.size == 1) {
          var out_domain = {0..out_shape[0]};
          var out_array: [out_domain] real = 0.0;
          out_dict.add(name, new DataArray1(out_array, property.dims));
        } else if (out_shape.size == 2) {
          var out_domain = {0..out_shape[0], 0..out_shape[1]};
          var out_array: [out_domain] real = 0.0;
          out_dict.add(name, new DataArray2(out_array, property.dims));
        } else if (out_shape.size == 3) {
          var out_domain = {0..out_shape[0], 0..out_shape[1], 0..out_shape[2]};
          var out_array: [out_domain] real = 0.0;
          out_dict.add(name, new DataArray3(out_array, property.dims));
        }
      }
    }

    return out_dict;
  }
  // var return_state = get_hybrid_sigma_pressure_levels();

  // var input_properties = new map(string, shared Properties);
  // input_properties.add("atmosphere_hybrid_sigma_pressure_a_coordinate_on_interface_levels", new shared Properties({"interface_levels", "*"}, new UnitMarker(0, 0, 0, 0, 0, 0, 0, 1, 0, "")));
  // input_properties.add("atmosphere_hybrid_sigma_pressure_b_coordinate_on_interface_levels", new shared Properties({"interface_levels", "*"}, new UnitMarker(0, 0, 0, 0, 0, 0, 0, 1, 0, "")));
  // input_properties.add("surface_air_pressure", new shared Properties({"*"}, new UnitMarker(-1, 1, -2, 0, 1, 0, 0, 1, 0, "Pa")));

  // var surface_arr: [0..27, 0..0] real = 1.0 * 1.0132e5;
  // return_state.add("surface_air_pressure", new shared DataArray2(surface_arr, {"lon", "lat"}));

  // var component = new HybridSigmaPressureDiagnosticComponent();
  // component.array_call(return_state);
  var slab = new shared SlabSurface();
  var components = [slab];
  var state = get_default_state(components);

  writeln(try! state.getValue("upwelling_shortwave_flux_in_air"));
}
