module RRTMGLongwave {
    use Components;
    use Map;
    use Properties;
    use marker;

    class RRTMGLongwave: TendencyComponent {
        proc init() {
            var input_properties = new map(string, shared Properties);
            input_properties.add("air_pressure", new shared Properties({"mid_levels"}, new UnitMarker(1, 0, 0, 0, 0, 0, 0, 1, 0, "m")));

            super.init(input_properties, new map(string, shared Properties), new map(string, shared Properties));
        }
    }
}