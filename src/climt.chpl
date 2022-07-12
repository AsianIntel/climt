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
  use Stepper;

  writeln("New library: climt");

  type dimensionless = UnitMarker(0, 0, 0, 0, 0, 0, 0);
  type pressure = UnitMarker(-1, 1, -2, 0, 0, 0, 0);
  type meter = UnitMarker(1, 0, 0, 0, 0, 0, 0);

  proc get_atmospheric_grid(name: string, value: real, state: State, marker: UnitMarker, interface_: bool = false, horizontal: bool = false) {
    var latitude = try! state.getValue("latitude"): DataArray2(real, dimensionless);
    var shape = latitude.arr.shape;

    var hybrid_sigma_pressure = try! state.getValue("atmosphere_hybrid_sigma_pressure_a_coordinate_on_interface_levels"): DataArray2(real, dimensionless);

    if horizontal {
      var arr: [0..shape[1]-1, 0..shape[0]-1] real = value;
      state.add(name, new shared DataArray2(arr, latitude.dimensions, marker));
    } else if interface_ {
      var arr: [0..hybrid_sigma_pressure.arr.shape[0]-1, 0..shape[1]-1, 0..shape[0]-1] real = value;
      var hsd = latitude.dimensions.these();
      var dom = {"interface_levels", hsd[0], hsd[1]};
      state.add(name, new shared DataArray3(arr, dom, marker));
    } else {
      var arr: [0..hybrid_sigma_pressure.arr.shape[0]-2, 0..shape[1]-1, 0..shape[0]-1] real = value;
      var hsd = latitude.dimensions.these();
      var dom = {"mid_levels", hsd[0], hsd[1]};
      state.add(name, new shared DataArray3(arr, dom, marker));
    }
  }

  proc get_surface_grid(name: string, value, state: State, marker: UnitMarker) {
    var latitude = try! state.getValue("latitude"): DataArray2(real, dimensionless);
    var shape = latitude.arr.shape;

    var arr: [0..shape[1]-1, 0..shape[0]-1] value.type = value;
    state.add(name, new shared DataArray2(arr, latitude.dimensions, marker));
  }

  proc get_land_grid(name: string, value: real, state: State, marker: UnitMarker, horizontal: bool = false) {
    var latitude = try! state.getValue("latitude"): DataArray2(real, dimensionless);
    var shape = latitude.arr.shape;

    if horizontal {
      var arr: [0..shape[1]-1, 0..shape[0]-1] value.type = value;
      state.add(name, new shared DataArray2(arr, latitude.dimensions, marker));
    }
  }

  proc get_ocean_grid(name: string, value: real, state: State, marker: UnitMarker, horizontal: bool = false) {
    var latitude = try! state.getValue("latitude"): DataArray2(real, dimensionless);
    var shape = latitude.arr.shape;

    if horizontal {
      var arr: [0..shape[1]-1, 0..shape[0]-1] real = value;
      state.add(name, new shared DataArray2(arr, latitude.dimensions, marker));
    }
  }

  proc get_ice_grid(name: string, value: real, state: State, marker: UnitMarker, interface_: bool = false, horizontal: bool = false) {
    var latitude = try! state.getValue("latitude"): DataArray2(real, dimensionless);
    var shape = latitude.arr.shape;

    var height_on_ice = try! state.getValue("height_on_ice_interface_levels"): DataArray1(real, meter);

    if horizontal {
      var arr: [0..shape[1]-1, 0..shape[0]-1] real = value;
      state.add(name, new shared DataArray2(arr, latitude.dimensions, marker));
    } else if interface_ {
      var arr: [0..height_on_ice.arr.shape[0]-1, 0..shape[1]-1, 0..shape[0]-1] real = value;
      var hsd = latitude.dimensions.these();
      var dom = {"ice_interface_levels", hsd[0], hsd[1]};
      state.add(name, new shared DataArray3(arr, dom, marker));
    } else {
      var arr: [0..height_on_ice.arr.shape[0]-2, 0..shape[1]-1, 0..shape[0]-1] real = value;
      var hsd = latitude.dimensions.these();
      var dom = {"ice_mid_levels", hsd[0], hsd[1]};
      state.add(name, new shared DataArray3(arr, dom, marker));
    }
  }

  proc default_values(name: string, state: State) {
    type velocity = UnitMarker(1, 0, -1, 0, 0, 0, 0);
    type time_inverse = UnitMarker(0, 0, -1, 0, 0, 0, 0);
    type mass_content = UnitMarker(-2, 1, 0, 0, 0, 0, 0);
    type temperature = UnitMarker(0, 0, 1, 0, 0, 0, 0);
    type length = UnitMarker(1, 0, 0, 0, 0, 0, 0);
    type flux = UnitMarker(0, 1, -3, 0, 0, 0, 0);

    select name {
      when "air_temperature" {
        get_atmospheric_grid(name, 290.0, state, new temperature(1, 0, "degK"), false, false);
      } 
      when "northward_wind" {
        get_atmospheric_grid(name, 0.0, state, new velocity(1, 0, "m/s"), false, false);
      } 
      when "eastward_wind" {
        get_atmospheric_grid(name, 0.0, state, new velocity(1, 0, "m/s"), false, false);
      } 
      when "divergence_of_wind" {
        get_atmospheric_grid(name, 0.0, state, new time_inverse(1, 0, "s-1"), false, false);
      } 
      when "atmospheric_relative_vorticity" {
        get_atmospheric_grid(name, 0.0, state, new time_inverse(1, 0, "s-1"), false, false);
      } 
      when "specific_humidity" {
        get_atmospheric_grid(name, 0.0, state, new dimensionless(1, 0, "kg/kg"), false, false);
      } 
      when "mole_fraction_of_carbon_dioxide_in_air" {
        get_atmospheric_grid(name, 330e-6, state, new dimensionless(1, 0, ""), false, false);
      } 
      when "mole_fraction_of_methane_in_air" {
        get_atmospheric_grid(name, 0.0, state, new dimensionless(1, 0, ""), false, false);
      } 
      when "mole_fraction_of_nitrous_oxide_in_air" {
        get_atmospheric_grid(name, 0.0, state, new dimensionless(1, 0, ""), false, false);
      } 
      when "mole_fraction_of_oxygen_in_air" {
        get_atmospheric_grid(name, 0.21, state, new dimensionless(1, 0, ""), false, false);
      } 
      when "mole_fraction_of_nitrogen_in_air" {
        get_atmospheric_grid(name, 0.78, state,  new dimensionless(1, 0, ""), false, false);
      } 
      when "mole_fraction_of_hydrogen_in_air" {
        get_atmospheric_grid(name, 500e-9, state,  new dimensionless(1, 0, ""), false, false);
      } 
      when "mole_fraction_of_cfc11_in_air" {
        get_atmospheric_grid(name, 0.0, state, new dimensionless(1, 0, ""), false, false);
      } 
      when "mole_fraction_of_cfc12_in_air" {
        get_atmospheric_grid(name, 0.0, state, new dimensionless(1, 0, ""), false, false);
      } 
      when "mole_fraction_of_cfc22_in_air" {
        get_atmospheric_grid(name, 0.0, state, new dimensionless(1, 0, ""), false, false);
      } 
      when "mole_fraction_of_carbon_tetrachloride_in_air" {
        get_atmospheric_grid(name, 0.0, state, new dimensionless(1, 0, ""), false, false);
      } 
      when "cloud_area_fraction_in_atmosphere_layer" {
        get_atmospheric_grid(name, 0.0, state, new dimensionless(1, 0, ""), false, false);
      } 
      when "mass_content_of_cloud_ice_in_atmosphere_layer" {
        get_atmospheric_grid(name, 0.0, state, new mass_content(1, 0, ""), false, false);
      } 
      when "mass_content_of_cloud_liquid_water_in_atmosphere_layer" {
        get_atmospheric_grid(name, 0.0, state, new mass_content(1, 0, ""), false, false);
      } 
      when "cloud_ice_particle_size" {
        get_atmospheric_grid(name, 20.0, state, new length(1000000, 0, ""), false, false);
      } 
      when "cloud_water_droplet_radius" {
        get_atmospheric_grid(name, 10.0, state, new length(1000000, 0, ""), false, false);
      } 
      when "cloud_base_mass_flux" {
        get_atmospheric_grid(name, 0.0, state, new UnitMarker(-2, 1, -1, 0, 0, 0, 0, 1, 0, "kg m^-2 s^-1"), false, true);
      } 
      when "zenith_angle" {
        get_atmospheric_grid(name, 0.0, state, new dimensionless(1, 0, "radians"), false, true);
      } 
      when "downwelling_shortwave_flux_in_air" {
        get_atmospheric_grid(name, 0.0, state, new flux(1, 0, "W m^-2"), true, false);
      } 
      when "downwelling_longwave_flux_in_air" {
        get_atmospheric_grid(name, 0.0, state, new flux(1, 0, "W m^-2"), true, false);
      } 
      when "upwelling_shortwave_flux_in_air" {
        get_atmospheric_grid(name, 0.0, state, new flux(1, 0, "W m^-2"), true, false);
      } 
      when "upwelling_longwave_flux_in_air" {
        get_atmospheric_grid(name, 0.0, state, new flux(1, 0, "W m^-2"), true, false);
      } 
      when "surface_specific_humidity" {
        get_surface_grid(name, 0.0, state, new dimensionless(1, 0, ""));
      } 
      when "surface_temperature" {
        get_surface_grid(name, 300.0, state, new temperature(1, 0, "degK"));
      } 
      when "soil_surface_temperature" {
        get_surface_grid(name, 300.0, state, new temperature(1, 0, "degK"));
      } 
      when "surface_geopotential" {
        get_surface_grid(name, 0.0, state, new UnitMarker(2, 0, -2, 0, 0, 0, 0, 1, 0, "m^2 s^-2"));
      } 
      when "surface_thermal_capacity" {
        get_surface_grid(name, 4.1813e3, state, new UnitMarker(2, 0, -2, 0, -1, 0, 0, 1, 0, "J kg^-1 degK^-1"));
      } 
      when "depth_of_slab_surface" {
        get_surface_grid(name, 50.0, state, new length(1, 0, "m"));
      } 
      when "surface_material_density" {
        get_surface_grid(name, 1000.0, state, new UnitMarker(-3, 1, 0, 0, 0, 0, 0, 1, 0, "kg m^-3"));
      } 
      when "surface_albedo_for_direct_shortwave" {
        get_surface_grid(name, 0.06, state, new dimensionless(1, 0, ""));
      } 
      when "surface_albedo_for_diffuse_shortwave" {
        get_surface_grid(name, 0.06, state, new dimensionless(1, 0, ""));
      } 
      when "surface_albedo_for_direct_near_infrared" {
        get_surface_grid(name, 0.06, state, new dimensionless(1, 0, ""));
      } 
      when "surface_albedo_for_diffuse_near_infrared" {
        get_surface_grid(name, 0.06, state, new dimensionless(1, 0, ""));
      } 
      when "surface_roughness_length" {
        get_surface_grid(name, 0.0002, state, new dimensionless(1, 0, ""));
      } 
      when "surface_drag_coefficient_for_heat_in_air" {
        get_surface_grid(name, 0.0012, state, new dimensionless(1, 0, ""));
      } 
      when "surface_drag_coefficient_for_momentum_in_air" {
        get_surface_grid(name, 0.0012, state, new dimensionless(1, 0, ""));
      } 
      when "area_type" {
        get_surface_grid(name, 0: int, state, new dimensionless(1, 0, ""));
      } 
      when "surface_upward_sensible_heat_flux" {
        get_surface_grid(name, 0.0, state, new flux(1, 0, "W m^-2"));
      } 
      when"surface_upward_latent_heat_flux" {
        get_surface_grid(name, 0.0, state, new flux(1, 0, "W m^-2"));
      } 
      when "soil_type" {
        get_land_grid(name, 0: int, state, new dimensionless(1, 0, ""), true);
      } 
      when "soil_temperature" {
        get_land_grid(name, 274.0, state, new temperature(1, 0, "degK"), true);
      } 
      when "soil_layer_thickness" {
        get_land_grid(name, 50.0, state, new length(1, 0, "m"), true);
      } 
      when "upward_heat_flux_at_ground_level_in_soil" {
        get_land_grid(name, 0.0, state, new flux(1, 0, "W m^-2"), true);
      } 
      when "heat_capacity_of_soil" {
        get_land_grid(name, 2000.0, state, new UnitMarker(2, 0, -2, 0, -1, 0, 0, 1, 0, "J kg^-1 degK^-1"), true);
      } 
      when "sea_water_density" {
        get_ocean_grid(name, 1.029e3, state, new UnitMarker(-3, 1, 0, 0, 0, 0, 0, 1, 0, "kg m^-3"), true);
      } 
      when "sea_surface_temperature" {
        get_ocean_grid(name, 300.0, state, new temperature(1, 0, "degK"), true);
      } 
      when "ocean_mixed_layer_thickness" {
        get_ocean_grid(name, 50.0, state, new length(1, 0, "m"), true);
      } 
      when "snow_and_ice_temperature" {
        get_ice_grid(name, 270.0, state, new temperature(1, 0, "degK"), true, false);
      } 
      when "heat_flux_into_sea_water_due_to_sea_ice" {
        get_ice_grid(name, 0.0, state, new flux(1, 0, "W m^-2"), false, true);
      }
      when "land_ice_thickness" {
        get_ice_grid(name, 0.0, state, new length(1, 0, "m"), false, true);
      } 
      when "sea_ice_thickness" {
        get_ice_grid(name, 0.0, state, new length(1, 0, "m"), false, true);
      } 
      when "surface_snow_thickness" {
        get_ice_grid(name, 0.0, state, new length(1, 0, "m"), false, true);
      } 
      when "solar_cycle_fraction" {
        var arr = [0.0];
        var dom: domain(string);
        state.add(name, new shared DataArray1(arr, dom, new dimensionless(1, 0, "")));
      } 
      when "flux_adjustment_for_earth_sun_distance" {
        var arr = [1.0];
        var dom: domain(string);
        state.add(name, new shared DataArray1(arr, dom, new dimensionless(1, 0, "")));
      } 
      when "lwe_thickness_of_soil_moisture_content" {
        get_surface_grid(name, 0.0, state, new length(1, 0, "m"));
      } 
      when "convective_precipitation_rate" {
        get_surface_grid(name, 0.0, state, new UnitMarker(1, 0, -1, 0, 0, 0, 0, 86400000, 0, "mm day^-1"));
      } 
      when "stratiform_precipitation_rate" {
        get_surface_grid(name, 0.0, state, new UnitMarker(1, 0, -1, 0, 0, 0, 0, 1, 0, "m s^-1"));
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
    return_state.add("surface_air_pressure", new shared DataArray2(surface_arr, {y_name, x_name}, new pressure(1, 0, "Pa")));

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
    return_state.add("longitude", new shared DataArray2(two_dim_lons, {y_name, x_name}, new dimensionless(1, 0, "")));

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
      return_state.add("latitude", new shared DataArray2(two_dim_lats, {y_name, x_name}, new dimensionless(1, 0, "")));
    } else if latitude_grid == "gaussian" {
      // TODO
    }

    var ice_interface_levels: [0..n_ice_interface_levels-1] real = 0.0;
    return_state.add("height_on_ice_interface_levels", new shared DataArray1(ice_interface_levels, {"ice_interface_levels"}, new meter(1, 0, "m")));

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
    state.add("atmosphere_hybrid_sigma_pressure_a_coordinate_on_interface_levels", new shared DataArray2(ak2, {"interface_levels", "*"}, new dimensionless(1, 0, "")));
    state.add("atmosphere_hybrid_sigma_pressure_b_coordinate_on_interface_levels", new shared DataArray2(bk2, {"interface_levels", "*"}, new dimensionless(1, 0, "")));

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

    var properties = new map(string, shared AbstractProperties);
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

  proc get_dim_lengths_from_raw_input(raw_input: State, input_properties: map(string, shared AbstractProperties)) throws {
    var dim_lengths = new map(string, int);

    for (name, property) in input_properties.items() {
        var i = 0;
        var ada = raw_input.getValue(name);
        for dim_name in property.getDims() {
          dim_lengths.add(dim_name, ada.getDomShape()[i]);
          i += 1;
        }   
    }

    return dim_lengths;
  }

  proc initialize_arrays_with_properties(output_properties: map(string, shared AbstractProperties), input_state: State, input_properties: map(string, shared AbstractProperties)) {
    var dim_lengths = get_dim_lengths_from_raw_input(input_state, input_properties);

    var out_dict = new map(string, shared AbstractDataArray);
    for (name, property) in output_properties.items() {
      var out_shape: [property.getDims().size] int;
      var i = 0;
      var marker = property.getMarker();
      for dim in property.getDims(){
        out_shape[i] = dim_lengths[dim];
        if (out_shape.size == 1) {
          var out_domain = {0..out_shape[0]};
          var out_array: [out_domain] real = 0.0;
          out_dict.add(name, new shared DataArray1(out_array, property.getDims(), marker));
        } else if (out_shape.size == 2) {
          // var out_domain = {0..out_shape[0], 0..out_shape[1]};
          // var out_array: [out_domain] real = 0.0;
          // out_dict.add(name, new shared DataArray2(out_array, property.getDims(), marker));
        } else if (out_shape.size == 3) {
          // var out_domain = {0..out_shape[0], 0..out_shape[1], 0..out_shape[2]};
          // var out_array: [out_domain] real = 0.0;
          // out_dict.add(name, new shared DataArray3(out_array, property.getDims(), marker));
        }
        i += 1;
      }
    }

    return out_dict;
  }

  var slab = new shared SlabSurface();
  var components = [slab];
  var state = get_default_state(components);
  var timestep = 3*60*60;

  var time_stepper = new AdamsBashforth(components);
  var (diagnostics, new_state) = time_stepper.call(state, timestep);
}
