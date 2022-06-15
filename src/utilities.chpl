module utilities {
    import Math;

    proc linspace(start: real, stop: real, num: int = 50, endpoint: bool = true) {
        var div;
        if (endpoint) {
            div = num - 1;
        } else {
            div = num;
        }
        var interval = (stop - start) / div;
        var array: [0..num-1] real;
        var val = start;
        for i in 0..num-1 {
            array[i] = val;
            val += interval;
        }

        return array;
    }

    proc sin(in arr) {
        for i in arr.domain {
            arr[i] = Math.sin(arr[i]);
        }

        return arr;
    }

    proc sum(in arr) {
        var sum: arr.eltType = 0;
        forall x in arr with (+ reduce sum) {
            sum += x;
        }
        return sum;
    }
}