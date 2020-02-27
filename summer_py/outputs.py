from matplotlib import pyplot
from matplotlib.ticker import FuncFormatter
import numpy
import os
import copy

from . import summer_model as sm
from . import post_processing as post_proc


def find_subplot_grid(n_plots):
    """
    Find a convenient number of rows and columns for a required number of subplots. First take the root of the number of
    subplots and round up to find the smallest square that could accommodate all of them. Next find out how many rows
    that many subplots would fill out by dividing the number of plots by the number of columns and rounding up. This
    will potentially leave a few panels blank at the end and number of rows will equal the number of columns or the
    number of rows will be on fewer.

    Args:
        n_plots: The number of subplots needed
    Returns:
        The number of rows of subplots
        n_cols: The number of columns of subplots
    """

    n_cols = int(numpy.ceil(numpy.sqrt(n_plots)))
    return int(numpy.ceil(n_plots / float(n_cols))), n_cols


def find_panel_grid_indices(axes, index, n_rows, n_columns):
    """
    Find the subplot index for a plot panel from the number of the panel and the number of columns of sub-plots.

    Args:
        axes: All the plot axes to be searched from
        index: The number of the panel counting up from zero
        n_rows: Number of rows of sub-plots in figure
        n_columns: Number of columns of sub-plots in figure
    """

    row, column = (
        numpy.floor_divide(index, n_columns),
        (index + 1) % n_columns - 1 if n_rows > 1 else None,
    )
    return axes[row, column] if n_rows > 1 else axes[index]


def get_label_font_size(max_dim):
    """
    Find standardised font size that can be applied across all figures.

    Args:
        max_dim: The number of rows or columns, whichever is the greater
    """

    label_font_sizes = {1: 8, 2: 7}
    return label_font_sizes[max_dim] if max_dim in label_font_sizes else 6


def scale_axes(vals, max_val, y_sig_figs):
    """
    General function to scale a set of axes and produce text that can be added to the axis label. Written here as a
    separate function from the tidy_axis method below because it can then be applied to both x- and y-axes.

    Args:
        vals: List of the current y-ticks
        max_val: The maximum value of this list
        y_sig_figs: The preferred number of significant figures for the ticks
    Returns:
        labels: List of the modified tick labels
        axis_modifier: The text to be added to the axis
    """

    y_number_format = "%." + str(y_sig_figs) + "f"
    y_number_format_around_one = "%." + str(max(2, y_sig_figs)) + "f"
    if max_val < 5e-9:
        labels = [y_number_format % (v * 1e12) for v in vals]
        axis_modifier = "Trillionth "
    elif max_val < 5e-6:
        labels = [y_number_format % (v * 1e9) for v in vals]
        axis_modifier = "Billionth "
    elif max_val < 5e-3:
        labels = [y_number_format % (v * 1e6) for v in vals]
        axis_modifier = "Millionth "
    elif max_val < 5e-2:
        labels = [y_number_format % (v * 1e3) for v in vals]
        axis_modifier = "Thousandth "
    elif max_val < 0.1:
        labels = [y_number_format % (v * 1e2) for v in vals]
        axis_modifier = "Hundredth "
    elif max_val < 5:
        labels = [y_number_format_around_one % v for v in vals]
        axis_modifier = ""
    elif max_val < 5e3:
        labels = [y_number_format % v for v in vals]
        axis_modifier = ""
    elif max_val < 5e6:
        labels = [y_number_format % (v / 1e3) for v in vals]
        axis_modifier = "Thousand "
    elif max_val < 5e9:
        labels = [y_number_format % (v / 1e6) for v in vals]
        axis_modifier = "Million "
    else:
        labels = [y_number_format % (v / 1e9) for v in vals]
        axis_modifier = "Billion "
    return labels, axis_modifier


def initialise_figures_axes(
    n_panels, room_for_legend=False, requested_grid=None, share_yaxis="none"
):
    """
    Initialise the subplots (or single plot) according to the number of panels required.

    Args:
        n_panels: The number of panels needed
        room_for_legend: Whether room is needed for a legend - applies to single axis plots only
        requested_grid: Shape of grid panels requested at call to method
        share_yaxis: String to pass to the sharey option
    Returns:
        fig: The figure object
        axes: A list containing each of the axes
        max_dims: The number of rows or columns of sub-plots, whichever is greater
    """

    pyplot.style.use("ggplot")
    n_rows, n_cols = requested_grid if requested_grid else find_subplot_grid(n_panels)
    horizontal_position_one_axis = 0.11 if room_for_legend else 0.15
    if n_panels == 1:
        fig = pyplot.figure()
        axes = fig.add_axes([horizontal_position_one_axis, 0.15, 0.69, 0.7])
    elif n_panels == 2:
        fig, axes = pyplot.subplots(1, 2)
        fig.set_figheight(3.5)
        fig.subplots_adjust(bottom=0.15, top=0.85)
    else:
        fig, axes = pyplot.subplots(n_rows, n_cols, sharey=share_yaxis)
        for panel in range(n_panels, n_rows * n_cols):
            find_panel_grid_indices(axes, panel, n_rows, n_cols).axis("off")
    return fig, axes, max([n_rows, n_cols]), n_rows, n_cols


def increment_list_for_patch(new_data, cumulative_data):
    """
    Takes a list of cumulative data totals, preserves the previous values and adds a new list to it. This is to allow
    patches to be plotted that have the previous data values as their base and the results of this stacking as their
    top.

    Args:
        new_data: The new data to be stacked up
        cumulative_data: The previous running totals
    Returns:
        previous_data: The previous running total (was cumulative_data)
        The new running total as the new values for cumulative_data
    """

    previous_data = copy.copy(cumulative_data)
    return previous_data, [last + current for last, current in zip(cumulative_data, new_data)]


def add_title_to_plot(fig, n_panels, content):
    """
    Function to add title to the top of a figure and handle multiple panels if necessary.

    Args:
        fig: The figure object to have a title added to it
        n_panels: Integer for the total number of panels on the figure
        content: Unprocessed string to determine text for the title
    """

    # if few panels, bigger and lower title
    greater_heights = {1: 0.92, 2: 0.98}
    greater_font_sizes = {1: 14, 2: 11}
    fig.suptitle(
        content,
        y=greater_heights[n_panels] if n_panels in greater_heights else 0.96,
        fontsize=greater_font_sizes[n_panels] if n_panels in greater_font_sizes else 10,
    )


class Outputs:
    def __init__(
        self,
        post_processing_list,
        targets_to_plot={},
        out_dir="outputs",
        translation_dict={},
        multiplot_only=False,
        mcmc_weights=None,
        plot_start_time=1990,
    ):
        """
        :param post_processing: an object of class post_processing associated with a run model
        :param out_dir: the name of the directory where to write the outputs
        """
        self.post_processing_list = post_processing_list
        self.targets_to_plot = targets_to_plot
        self.out_dir = out_dir
        self.translation_dict = translation_dict
        self.multiplot_only = multiplot_only
        self.scenario_names = {}
        self.plot_start_time = plot_start_time

        self.colour_theme = [
            (0.0, 0.0, 0.0),
            (57.0 / 255.0, 106.0 / 255.0, 177.0 / 255.0),
            (218.0 / 255.0, 124.0 / 255.0, 48.0 / 255.0),
            (62.0 / 255.0, 150.0 / 255.0, 81.0 / 255.0),
            (204.0 / 255.0, 37.0 / 255.0, 41.0 / 255.0),
            (107.0 / 255.0, 76.0 / 255.0, 154.0 / 255.0),
            (146.0 / 255.0, 36.0 / 255.0, 40.0 / 255.0),
            (148.0 / 255.0, 139.0 / 255.0, 61.0 / 255.0),
            (0.0, 0.0, 125.0 / 255.0),
            (210.0 / 255.0, 70.0 / 255.0, 0.0),
            (100.0 / 255.0, 150.0 / 255.0, 1.0),
            (65.0 / 255.0, 65.0 / 255.0, 65.0 / 255.0),
            (220.0 / 255.0, 25.0 / 255.0, 25.0 / 255.0),
            (120.0 / 255.0, 55.0 / 255.0, 20.0 / 255.0),
            (120.0 / 255.0, 55.0 / 255.0, 110.0 / 255.0),
            (135.0 / 255.0, 135.0 / 255.0, 30.0 / 255.0),
            (120.0 / 255.0, 120.0 / 255.0, 120.0 / 255.0),
            (220.0 / 255.0, 20.0 / 255.0, 170.0 / 255.0),
            (20.0 / 255.0, 65.0 / 255.0, 20.0 / 255.0),
            (15.0 / 255.0, 145.0 / 255.0, 25.0 / 255.0),
            (15.0 / 255.0, 185.0 / 255.0, 240.0 / 255.0),
            (10.0 / 255.0, 0.0, 110.0 / 255.0),
            (0.5, 0.5, 0.5),
            (0.0, 0.0, 0.0),
        ]

        self.mcmc_weights = mcmc_weights
        self.mcmc_mode = not self.mcmc_weights is None
        self.multiplot_only = True is self.mcmc_mode

        self.create_out_directories()

    def create_out_directories(self):
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)

        for scenario_index in range(len(self.post_processing_list)):
            scenario_number = self.post_processing_list[scenario_index].scenario_number
            scenario_name = (
                "Baseline" if scenario_number == 0 else "Scenario " + str(scenario_number)
            )
            self.scenario_names[scenario_number] = scenario_name

            if not self.multiplot_only and not self.mcmc_mode:
                scenario_out_dir = os.path.join(self.out_dir, scenario_name)
                if not os.path.exists(scenario_out_dir):
                    os.mkdir(scenario_out_dir)

        if len(self.scenario_names) > 1 or self.mcmc_mode:
            multi_out_dir = os.path.join(self.out_dir, "multi_plots")
            if not os.path.exists(multi_out_dir):
                os.mkdir(multi_out_dir)

    def tidy_x_axis(self, axis, start, end, max_dims, labels_off=False, x_label=None):
        """
        Function to tidy x-axis of a plot panel - currently only used in the scale-up vars, but intended to be written
        in such a way as to be extendable to other types of plotting.

        Args:
            axis: The plotting axis
            start: Lowest x-value being plotted
            end: Highest x-value being plotted
            max_dim: Maximum number of rows or columns of subplots in figure
            labels_off: Whether to turn all tick labels off on this axis
            x_label: Text for the x-axis label if required
        """

        # range
        axis.set_xlim(left=start, right=end)

        # ticks and their labels
        if labels_off:
            axis.tick_params(axis="x", labelbottom=False)
        elif len(axis.get_xticks()) > 7:
            for label in axis.xaxis.get_ticklabels()[::2]:
                label.set_visible(False)
        axis.tick_params(axis="x", length=3, pad=6, labelsize=get_label_font_size(max_dims))

        # axis label
        if x_label is not None:
            axis.set_xlabel(
                self.intelligent_convert_string(x_label), fontsize=get_label_font_size(max_dims)
            )

    def tidy_y_axis(
        self,
        axis,
        quantity,
        max_dims,
        left_axis=True,
        max_value=1e6,
        space_at_top=0.1,
        y_label=None,
        y_lims=None,
        allow_negative=False,
    ):
        """
        General approach to tidying up the vertical axis of a plot, depends on whether it is the left-most panel.

        Args:
            axis: The axis itself
            quantity: The name of the quantity being plotted (which can be used to determine the sort of variable it is)
            max_dims: Maximum number of rows or columns of subplots on the figure
            left_axis: Boolean for whether the axis is the left-most panel
            max_value: The maximum value in the data being plotted
            space_at_top: Relative amount of space to leave at the top, above the maximum value of the plotted data
            y_label: A label for the y-axis, if required
            y_lims: 2-element tuple for the y-limit, if required
            allow_negative: Whether to set the bottom of the axis to zero
        """

        # axis range
        if y_lims:
            axis.set_ylim(y_lims)
        elif "prop_" in quantity and axis.get_ylim()[1] > 1.0:
            axis.set_ylim(top=1.004)
        elif "prop_" in quantity and max_value > 0.7:
            axis.set_ylim(bottom=0.0, top=1.0)
        elif "prop_" in quantity or "likelihood" in quantity or "cost" in quantity:
            pass
        elif axis.get_ylim()[1] < max_value * (1.0 + space_at_top):
            pass
            # axis.set_ylim(top=max_value * (1. + space_at_top))
        if not allow_negative:
            axis.set_ylim(bottom=0.0)

        # ticks
        axis.tick_params(axis="y", length=3.0, pad=6, labelsize=get_label_font_size(max_dims))

        # tick labels
        if not left_axis:
            pyplot.setp(axis.get_yticklabels(), visible=False)
        elif "prop_" in quantity:
            axis.yaxis.set_major_formatter(FuncFormatter("{0:.0%}".format))

        # axis label
        if y_label and left_axis:
            axis.set_ylabel(
                self.intelligent_convert_string(y_label), fontsize=get_label_font_size(max_dims)
            )

    def intelligent_convert_string(self, string):
        """
        returns a more readable version of string for figure title ...
        """
        if string in self.translation_dict.keys():
            return self.translation_dict[string]
        elif string[0:22] == "distribution_of_strata":
            return "Population distribution by " + string.split("X")[1]
        elif string[0:4] == "prev":
            char = "Prevalence of "
            numerator_groups = string.split("among")[0].split("X")[1:]
            for group in numerator_groups:
                char += group + " "

            subgroup = string.split("among")[1]
            if len(subgroup) > 0:
                char += "("
                need_comma = False
                for group in subgroup.split("X"):
                    if len(group) > 0:
                        if need_comma:
                            char += ", "
                        char += group
                        need_comma = True
                char += ")"
            return char
        else:
            return string

    def finish_off_figure(self, fig, filename, n_plots=1, title_text=None):
        """
        Slight extension of save_figure to include adding main title to figure.

        Args:
            fig: The figure to add the title to
            filename: The end of the string for the file name
            n_plots: number of panels
            title_text: Text for the title of the figure
        """
        if title_text is not None:
            add_title_to_plot(fig, n_plots, self.intelligent_convert_string(title_text))
        filename = os.path.join(self.out_dir, filename + ".png")
        fig.savefig(filename, dpi=300, bbox_inches="tight")
        pyplot.close(fig)

    def plot_requested_outputs(self):
        """
        main method to run the plotting of all the outputs requested in the post-processing object
        """
        outputs_to_plot = (
            self.post_processing_list[0].requested_outputs
            if self.post_processing_list[0].derived_outputs is None
            else self.post_processing_list[0].requested_outputs
            + list(self.post_processing_list[0].derived_outputs.keys())
        )

        for requested_output in outputs_to_plot:
            if requested_output == "times":
                continue
            if self.multiplot_only:
                multiplot_plotting_modes = [True]
            elif len(self.scenario_names) > 1:
                multiplot_plotting_modes = [False, True]
            else:
                multiplot_plotting_modes = [False]

            for multi_plot in multiplot_plotting_modes:
                if requested_output not in self.post_processing_list[0].derived_outputs:
                    if (
                        isinstance(
                            self.post_processing_list[0].generated_outputs[requested_output], dict
                        )
                        and requested_output[0:22] == "distribution_of_strata"
                        and multi_plot
                    ):
                        continue

                y_max = -1.0e9

                if self.mcmc_mode:
                    mcmc_array = None  # initialise

                for scenario_index, scenario_number in enumerate(self.scenario_names.keys()):
                    data_to_plot = (
                        self.post_processing_list[scenario_index].generated_outputs[
                            requested_output
                        ]
                        if requested_output
                        not in self.post_processing_list[scenario_index].derived_outputs
                        else self.post_processing_list[scenario_index].derived_outputs[
                            requested_output
                        ]
                    )

                    if self.mcmc_mode:
                        this_series_array = numpy.array([data_to_plot])
                        if mcmc_array is None:
                            mcmc_array = copy.copy(this_series_array)
                        else:
                            mcmc_array = numpy.concatenate((mcmc_array, this_series_array), axis=0)

                    if not multi_plot or scenario_index == 0:
                        fig, axes, max_dims, n_rows, n_cols = initialise_figures_axes(1)
                        axis = find_panel_grid_indices([axes], 0, n_rows, n_cols)
                        output_name = requested_output

                        # plot targets
                        if requested_output in self.targets_to_plot.keys():
                            self.plot_targets(requested_output, axis)

                    times_to_plot = (
                        self.post_processing_list[scenario_index].model.times
                        if requested_output
                        not in self.post_processing_list[scenario_index].requested_times.keys()
                        else self.post_processing_list[scenario_index].requested_times[
                            requested_output
                        ]
                    )

                    if isinstance(data_to_plot, list):

                        if not multi_plot and scenario_index > 0:
                            times_to_plot_0 = (
                                self.post_processing_list[0].model.times
                                if requested_output
                                not in self.post_processing_list[0].requested_times.keys()
                                else self.post_processing_list[0].requested_times[requested_output]
                            )

                            data_to_plot_0 = (
                                self.post_processing_list[0].generated_outputs[requested_output]
                                if requested_output
                                not in self.post_processing_list[0].derived_outputs
                                else self.post_processing_list[0].derived_outputs[requested_output]
                            )

                            axis.plot(
                                times_to_plot_0,
                                data_to_plot_0,
                                color=self.colour_theme[0],
                                label="Baseline",
                            )

                        this_label = (
                            "Scenario " + str(scenario_number)
                            if scenario_number > 0
                            else "Baseline"
                        )

                        sc_color = (
                            self.colour_theme[scenario_number]
                            if scenario_number < len(self.colour_theme) and not self.mcmc_mode
                            else "black"
                        )

                        axis.plot(times_to_plot, data_to_plot, color=sc_color, label=this_label)

                        if requested_output in self.post_processing_list[scenario_index].ymax:
                            y_max = self.post_processing_list[scenario_index].ymax[requested_output]
                        else:
                            y_max = max([y_max, max(data_to_plot)])
                        if scenario_index == 0 or not multi_plot:
                            self.tidy_x_axis(
                                axis,
                                start=self.plot_start_time,
                                end=max(times_to_plot),
                                max_dims=max_dims,
                                x_label="time",
                            )

                        if not multi_plot or scenario_index == len(self.scenario_names) - 1:
                            self.tidy_y_axis(
                                axis,
                                quantity="",
                                max_dims=max_dims,
                                y_label=output_name,
                                max_value=y_max,
                            )

                    elif (
                        isinstance(data_to_plot, dict)
                        and requested_output[0:22] == "distribution_of_strata"
                    ):

                        current_data = data_to_plot
                        self.plot_stacked_epi_outputs(
                            axis, times_to_plot, current_data, fraction=False
                        )

                    if not multi_plot or scenario_index == len(self.scenario_names) - 1:
                        if multi_plot:
                            dir_name = "multi_plots"
                            if not self.mcmc_mode:
                                axis.legend(bbox_to_anchor=(1.0, 1))
                        else:
                            dir_name = self.scenario_names[scenario_number]
                            if scenario_index > 0:
                                axis.legend(bbox_to_anchor=(1.0, 1))

                        file_name = os.path.join(dir_name, output_name)
                        self.finish_off_figure(fig, filename=file_name, title_text=output_name)

                if self.mcmc_mode:
                    percentiles = numpy.percentile(mcmc_array, [2.5, 50, 97.5], axis=0)
                    fig, axes, max_dims, n_rows, n_cols = initialise_figures_axes(1)
                    axis = find_panel_grid_indices([axes], 0, n_rows, n_cols)
                    output_name = requested_output
                    for perc_index in range(3):
                        axis.plot(
                            times_to_plot,
                            percentiles[perc_index],
                            color="black",
                            linestyle="-" if perc_index == 1 else "--",
                        )
                    # plot targets
                    if requested_output in self.targets_to_plot.keys():
                        self.plot_targets(requested_output, axis)

                    if requested_output in self.post_processing_list[scenario_index].ymax:
                        y_max = self.post_processing_list[scenario_index].ymax[requested_output]
                    else:
                        y_max = percentiles.max()
                    self.tidy_x_axis(
                        axis,
                        start=self.plot_start_time,
                        end=max(times_to_plot),
                        max_dims=max_dims,
                        x_label="time",
                    )
                    self.tidy_y_axis(
                        axis, quantity="", max_dims=max_dims, y_label=output_name, max_value=y_max
                    )
                    dir_name = "multi_plots"
                    file_name = os.path.join(dir_name, output_name + "_ci")
                    self.finish_off_figure(fig, filename=file_name, title_text=output_name)

    def plot_targets(self, requested_output, axis):

        targets = self.targets_to_plot[requested_output]
        for i, time in enumerate(targets["times"]):
            if len(targets["values"][i]) > 1:  # plot CI
                x_vals = [time, time]
                y_vals = targets["values"][i][1:]
                axis.plot(x_vals, y_vals, "m", linewidth=1, color="red")

            # plot point estimate
            marker_size = 30.0
            for colour in ["red", "white"]:
                axis.scatter(
                    time, targets["values"][i][0], marker="o", color=colour, s=marker_size,
                )
                marker_size -= 20.0

    def plot_stacked_epi_outputs(self, axis, times_to_plot, current_data, fraction=True):
        # plot patches and proxy by category
        cumulative_data = [0.0] * len(times_to_plot)
        current_data_as_array = numpy.array(list(current_data.values()))
        populations = current_data_as_array.sum(0)

        if fraction:
            current_data_as_array = numpy.divide(current_data_as_array, populations)

        for l, stratum in enumerate(current_data):
            if fraction:
                current_data[stratum] = current_data_as_array[l, :]
            previous_data, cumulative_data = increment_list_for_patch(
                current_data[stratum], cumulative_data
            )
            colour = self.colour_theme[l + 1]

            if len(previous_data) == len(cumulative_data) and len(times_to_plot) == len(
                previous_data
            ):
                axis.fill_between(
                    times_to_plot,
                    previous_data,
                    cumulative_data,
                    facecolor=colour,
                    edgecolor=colour,
                    alpha=0.8,
                )
            axis.plot([-1e2], [0.0], color=colour, label=stratum, linewidth=5.0)  # proxy

        axis.legend(bbox_to_anchor=(1.0, 1))

        if fraction:
            y_label = "proportion"
            y_max = 1.0
        else:
            y_label = "population size"
            y_max = max(populations)

        self.tidy_x_axis(
            axis, start=self.plot_start_time, end=max(times_to_plot), max_dims=1, x_label="time"
        )
        self.tidy_y_axis(axis, quantity="", max_dims=1, y_label=y_label, max_value=y_max)

    def plot_outputs_by_stratum(self, requested_output="prevXinfectious", sc_index=0):
        if not hasattr(self.post_processing_list[sc_index].model, "all_stratifications"):
            return
        all_groups = self.post_processing_list[sc_index].model.all_stratifications
        for stratification in all_groups.keys():
            fig, axes, max_dims, n_rows, n_cols = initialise_figures_axes(1)
            axis = find_panel_grid_indices([axes], 0, n_rows, n_cols)

            times_to_plot = (
                self.post_processing_list[sc_index].model.times
                if requested_output
                not in self.post_processing_list[sc_index].requested_times.keys()
                else self.post_processing_list[sc_index].requested_times[requested_output]
            )

            times_to_plot_baseline = (
                self.post_processing_list[0].model.times
                if requested_output
                not in self.post_processing_list[sc_index].requested_times.keys()
                else self.post_processing_list[0].requested_times[requested_output]
            )

            for i, stratum in enumerate(all_groups[stratification]):
                requested_output_for_stratum = (
                    requested_output + "XamongX" + stratification + "_" + stratum
                )
                _label = (
                    self.translation_dict[stratification + "_" + stratum]
                    if stratification + "_" + stratum in self.translation_dict.keys()
                    else stratification + "_" + stratum
                )

                # plot baseline run in addition to scenario run
                if sc_index > 0:
                    axis.plot(
                        times_to_plot_baseline,
                        self.post_processing_list[0].generated_outputs[
                            requested_output_for_stratum
                        ],
                        color=self.colour_theme[0],
                        label="baseline",
                    )

                axis.plot(
                    times_to_plot,
                    self.post_processing_list[sc_index].generated_outputs[
                        requested_output_for_stratum
                    ],
                    color=self.colour_theme[i + 1],
                    label=_label,
                )

            self.tidy_x_axis(
                axis,
                start=self.plot_start_time,
                end=max(times_to_plot),
                max_dims=max_dims,
                x_label="time",
            )
            self.tidy_y_axis(
                axis, quantity="", max_dims=max_dims, y_label=requested_output + "Xamong"
            )

            axis.legend(bbox_to_anchor=(1.0, 1))

            scenario_name = list(self.scenario_names.values())[sc_index]

            file_name = os.path.join(scenario_name, requested_output + "BY" + stratification)
            self.finish_off_figure(fig, filename=file_name, title_text=requested_output + "Xamong")

    def plot_input_function(self, input_function_name, input_function, sc_index=0):
        """
        Plot single simple plot of a function
        """

        times_to_plot = self.post_processing_list[sc_index].model.times
        fig, axes, max_dims, n_rows, n_cols = initialise_figures_axes(1)
        axis = find_panel_grid_indices([axes], 0, n_rows, n_cols)
        axis.plot(
            times_to_plot,
            list(map(input_function, times_to_plot)),
            color=self.colour_theme[sc_index])
        self.tidy_x_axis(
            axis,
            start=self.plot_start_time,
            end=max(times_to_plot),
            max_dims=max_dims,
            x_label="time",
        )
        self.tidy_y_axis(
            axis, quantity="", max_dims=max_dims
        )
        scenario_name = list(self.scenario_names.values())[sc_index]
        file_name = os.path.join(scenario_name, input_function_name)
        self.finish_off_figure(fig, filename=file_name, title_text=input_function_name)


if __name__ == "__main__":
    # build and run an example model
    sir_model = sm.StratifiedModel(
        numpy.linspace(0, 60 / 365, 61).tolist(),
        ["susceptible", "infectious", "recovered"],
        {"infectious": 0.001},
        {"beta": 400, "recovery": 365 / 13, "infect_death": 1},
        [
            {
                "type": "standard_flows",
                "parameter": "recovery",
                "origin": "infectious",
                "to": "recovered",
            },
            {
                "type": "infection_frequency",
                "parameter": "beta",
                "origin": "susceptible",
                "to": "infectious",
            },
            {"type": "compartment_death", "parameter": "infect_death", "origin": "infectious"},
        ],
        output_connections={"incidence": {"origin": "susceptible", "to": "infectious"}},
        verbose=False,
        integration_type="solve_ivp",
    )

    sir_model.stratify(
        "strain",
        ["sensitive", "resistant"],
        ["infectious"],
        requested_proportions={},
        verbose=False,
    )

    age_mixing = None
    sir_model.stratify(
        "age",
        [1, 10, 3],
        [],
        {},
        {"recovery": {"1": 0.5, "10": 0.8}},
        infectiousness_adjustments={"1": 0.8},
        mixing_matrix=age_mixing,
        verbose=False,
    )

    sir_model.run_model()

    # request some outputs
    req_outputs = [
        "prevXinfectiousXamongXage_10Xstrain_sensitive",
        # 'prevXinfectiousXresistantXamongXage_10Xstrain_sensitive',
        "distribution_of_strataXstrain",
        "distribution_of_strataXage",
    ]
    multipliers = {
        "prevXinfectiousXamongXage_10Xstrain_sensitive": 1.0e5,
        "prevXinfectiousXamong": 1.0e5,
    }
    pp = post_proc.PostProcessing(
        sir_model, req_outputs, multipliers=multipliers, scenario_number=0
    )
    pp2 = post_proc.PostProcessing(
        sir_model, req_outputs, multipliers=multipliers, scenario_number=1
    )

    targets_to_plot = {
        "prevXinfectiousXamongXage_10Xstrain_sensitive": {
            "times": [0.01, 0.02],
            "values": [[20000, 18000, 22000], [23000]],
        }
    }
    translation_dict = {
        "prevXinfectiousXamongXage_10Xstrain_sensitive": "Prevalence of TB among 10-30 years old"
    }

    # generate outputs
    outputs = Outputs([pp, pp2], targets_to_plot=targets_to_plot, translation_dict=translation_dict)
    outputs.plot_requested_outputs()
