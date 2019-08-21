
def create_function_of_function(outer_function, inner_function):
    """
    function that can itself return a function that sequentially apply two functions

    :param outer_function: function
        last function to be called
    :param inner_function: function
        first function to be called
    :return: function
        composite function that applies the inner and then the outer function, allowing the time parameter to be passed
            through if necessary
    """

    def function_to_return(time):
        return outer_function(inner_function(time))
    return function_to_return


def inner_function():
    def inner_funct(time):
        return 0.5
    return inner_funct


def multiplier_function(multiplier):
    def multiplier_funct(input):
        return multiplier * input
    return multiplier_funct


multiplier = multiplier_function(10.0)
inner = inner_function

def get_composite_function(inner_function, outer_function):
    print(outer_function(inner_function(2000.)))

print(get_composite_function(inner, multiplier))


