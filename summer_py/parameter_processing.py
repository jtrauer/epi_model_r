

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
