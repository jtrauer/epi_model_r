import numpy


def change_parameter_unit(parameter_dict, multiplier):
    """
    used to adapt the latency parameters from the earlier functions according to whether they are needed as by year
        rather than by day
    :param parameter_dict: dict
        dictionary whose values need to be adjusted
    :param multiplier: float
        multiplier
    :return: dict
        dictionary with values multiplied by the multiplier argument
    """
    return {param_key: param_value * multiplier for param_key, param_value in parameter_dict.items()}


def add_w_to_param_names(parameter_dict):
    """
    add a "W" string to the end of the parameter name to indicate that we should over-write up the chain
    :param parameter_dict: dict
        the dictionary before the adjustments
    :return: dict
        same dictionary but with the "W" string added to each of the keys
    """
    return {str(age_group) + "W": value for age_group, value in parameter_dict.items()}


def create_step_function_from_dict(input_dict):
    """
    create a step function out of dictionary with numeric keys and values, where the keys determine the values of the
    independent variable at which the steps between the output values occur
    :param input_dict: dict
        dictionary in standard format as described above
    :return: function
        the function constructed from input data
    """
    dict_keys = list(input_dict.keys())
    dict_keys.sort()
    dict_values = [input_dict[key] for key in dict_keys]

    def step_function(input_value):
        if input_value >= dict_keys[-1]:
            return dict_values[-1]
        else:
            for key in range(len(dict_keys)):
                if input_value < dict_keys[key + 1]:
                    return dict_values[key]

    return step_function


def sinusoidal_scaling_function(start_time, baseline_value, end_time, final_value):
    """
    in order to implement scale-up functions over time, use the cosine function to produce smooth scale-up functions
    from one point to another, returning the starting value before the starting point and the final value after the
    end point

    :param start_time: float
        starting value of the independent variable
    :param baseline_value: float
        starting value of the dependent variable
    :param end_time: float
        final value of the independent variable
    :param final_value: float
        final value of the dependent variable
    :return:
        function scaling from the starting value to the final value
    """
    def sinusoidal_function(time):
        if not isinstance(time, float):
            raise ValueError("value provided to scaling function not a float")
        elif start_time > end_time:
            raise ValueError("start time is later than end time")
        elif time < start_time:
            return baseline_value
        elif start_time <= time <= end_time:
            return baseline_value + \
                   (final_value - baseline_value) * \
                   (0.5 - 0.5 * numpy.cos((time - start_time) * numpy.pi / (end_time - start_time)))
        else:
            return final_value
    return sinusoidal_function


def logistic_scaling_function(parameter):
    """
    a specific sigmoidal form of function that scales up from zero to one around the point of parameter
        won't be useful in all situations and is specifically for age-specific infectiousness - should be the same as in
        Romain's BMC Medicine manuscript

    :param parameter: float
        the single parameter to the function
    :return: function
        the logistic function
    """
    return lambda x: 1.0 - 1.0 / (1.0 + numpy.exp(-(parameter - x)))
