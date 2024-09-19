import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors

class MultiPlotter:
    def __init__(self, nrows: int = 1, ncols: int = 1, figsize: tuple = (10, 8), suptitle: str = None):
        """
        Initializes the MultiPlotter class.
        
        Parameters:
        - nrows (int): Number of rows of subplots.
        - ncols (int): Number of columns of subplots.
        - figsize (tuple): Size of the entire figure (width, height).
        - suptitle (str): Overall title for the entire figure.
        """
        self.fig, self.axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
        self.axes = self.axes.flatten() if nrows * ncols > 1 else [self.axes]
        self.current_ax = 0
        if suptitle:
            self.fig.suptitle(suptitle, fontsize=20)

    def plot_hist(self, data: np.ndarray,
                  bins: np.ndarray,
                  weights: np.ndarray = None,
                  log_scale: bool = False, 
                  label: str = '', 
                  color: str = 'blue',
                  xlabel: str = None,
                  ylabel: str = None,
                  xlim: tuple = None,
                  ylim: tuple = None,
                  ticks: np.ndarray = None):
        """
        Plots a 1D histogram on the current axis.

        Parameters:
        - data (np.ndarray): Data to plot in the histogram.
        - bins (np.ndarray): Bins of the histogram.
        - log_scale (bool): Whether to use a logarithmic scale for the y-axis.
        - label (str): Label of the histogram.
        """
        ax = self.axes[self.current_ax]
        if weights is None:
            weights = np.ones(len(data))
        ax.hist(data, bins=bins, weights=weights, log=log_scale, histtype='step', label=label, color=color)
        self._apply_formatting(ax, xlabel, ylabel, xlim, ylim, ticks)

    def _apply_formatting(self, ax, xlabel=None, ylabel=None, xlim=None, ylim=None, ticks=None):
        """
        Apply common formatting to the axis.
        """
        if xlabel: ax.set_xlabel(xlabel, fontsize=20)
        if ylabel: ax.set_ylabel(ylabel, fontsize=20)
        if xlim: ax.set_xlim(xlim)
        if ylim: ax.set_ylim(ylim)
        if ticks is not None: 
            ax.set_xticks(ticks)
            ax.tick_params(axis='x', labelsize=12)

    def plot_hist2d(self, 
                    x: np.ndarray, 
                    y: np.ndarray, 
                    bins_x: np.ndarray, 
                    bins_y: np.ndarray, 
                    log_scale: bool = False, 
                    xlabel: str = None,
                    ylabel: str = None,
                    xlim: tuple = None,
                    ylim: tuple = None,
                    ticks: np.ndarray = None):
        """
        Plots a 2D histogram on the current axis.

        Parameters:
        - x (np.ndarray): Data for the x-axis.
        - y (np.ndarray): Data for the y-axis.
        - bins_x (np.ndarray): Number of bins for x.
        - bins_y (np.ndarray): Number of bins for y.
        - log_scale (bool): Whether to use a logarithmic scale for the color scale.
        - xlabel (str): Label for the x-axis.
        - ylabel (str): Label for the y-axis.
        - xlim (tuple): Limits for the x-axis (min, max).
        - ylim (tuple): Limits for the y-axis (min, max).
        - ticks (np.ndarray): Ticks for the x-axis.
        """
        ax = self.axes[self.current_ax]
        hist = ax.hist2d(x, y, bins=(bins_x, bins_y), norm=colors.LogNorm() if log_scale else None)
        self.fig.colorbar(hist[3], ax=ax)
        self._apply_formatting(ax, xlabel, ylabel, xlim, ylim, ticks)


    def set_labels(self, xlabel: str = None, ylabel: str = None, fontsize: int = 20):
        """
        Sets the x and y labels for the current axis.
        
        Parameters:
        - xlabel (str): Label for the x-axis.
        - ylabel (str): Label for the y-axis.
        """
        ax = self.axes[self.current_ax]
        if(xlabel): ax.set_xlabel(xlabel, fontsize=fontsize)
        if(ylabel): ax.set_ylabel(ylabel, fontsize=fontsize)

    def set_limits(self, xlim: tuple = None, ylim: tuple = None):
        """
        Sets the x and y limits for the current axis.
        
        Parameters:
        - xlim (tuple): Limits for the x-axis (min, max).
        - ylim (tuple): Limits for the y-axis (min, max).
        """
        ax = self.axes[self.current_ax]
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
    
    def set_xticks(self, ticks : np.ndarray = None):
        ax = self.axes[self.current_ax]
        if ticks:
            ax.set_xticks(ticks)

    def add_grid(self, which: str = 'both', linestyle: str = '--', color: str = 'gray'):
        """
        Adds grid lines to the current axis.
        
        Parameters:
        - which (str): Which grid lines to show ('both', 'x', 'y').
        - linestyle (str): Style of the grid lines.
        - color (str): Color of the grid lines.
        """
        ax = self.axes[self.current_ax]
        ax.grid(which=which, linestyle=linestyle, color=color)

    def add_legend(self, labels: list):
        """
        Adds a legend to the current axis.
        
        Parameters:
        - labels (list): List of labels for the legend.
        """
        ax = self.axes[self.current_ax]
        ax.legend(labels, fontsize=15)

    def customize_ticks(self, xticks: list = None, yticks: list = None):
        """
        Customizes the ticks on the x and y axes.
        
        Parameters:
        - xticks (list): Custom ticks for the x-axis.
        - yticks (list): Custom ticks for the y-axis.
        """
        ax = self.axes[self.current_ax]
        if xticks:
            ax.set_xticks(xticks)
        if yticks:
            ax.set_yticks(yticks)

    def next_plot(self, plot_legend: bool = False):
        """
        Moves to the next subplot axis.
        Plot legend if specified.
        """
        if(plot_legend): self.axes[self.current_ax].legend()
        self.current_ax += 1

    def show(self):
        """
        Displays the entire figure with all subplots.
        """
        plt.tight_layout()
        plt.show()

    def save_to_pdf(self, filename: str):
        """
        Saves the entire figure with all subplots to a PDF file.
        
        Parameters:
        - filename (str): The name of the PDF file to save the figure.
        """
        with PdfPages(filename) as pdf:
            pdf.savefig(self.fig)
            plt.close(self.fig)
    
    def save_multiple_figures_to_pdf(self, filename: str, plotters: list):
        """
        Saves multiple figures (from different MultiPlotter instances) to a single PDF file.
        
        Parameters:
        - filename (str): The name of the PDF file to save the figures.
        - plotters (list): A list of MultiPlotter instances whose figures will be saved.
        """
        with PdfPages(filename) as pdf:
            for plotter in plotters:
                pdf.savefig(plotter.fig)
                plt.close(plotter.fig)

    # Decorators for adding common features easily
    def with_title(title: str):
        """
        Decorator to add a title to the current subplot.
        
        Parameters:
        - title (str): Title to add to the subplot.
        """
        def decorator(func):
            def wrapper(self, *args, **kwargs):
                func(self, *args, **kwargs)
                ax = self.axes[self.current_ax]
                ax.set_title(title, fontsize=16)
            return wrapper
        return decorator

    def with_grid(linestyle: str = '--', color: str = 'gray'):
        """
        Decorator to add grid lines to the current subplot.
        
        Parameters:
        - linestyle (str): Style of the grid lines.
        - color (str): Color of the grid lines.
        """
        def decorator(func):
            def wrapper(self, *args, **kwargs):
                func(self, *args, **kwargs)
                self.add_grid(linestyle=linestyle, color=color)
            return wrapper
        return decorator

# # Example usage
# if __name__ == "__main__":
#     plot_tools = MultiPlotter(nrows=1, ncols=2, figsize=(14, 6), suptitle="1D and 2D Histograms with Decorators")

#     @plot_tools.with_title("1D Histogram")
#     @plot_tools.with_grid()
#     def plot_first_hist(plot_tools):
#         data1 = np.random.normal(size=1000)
#         plot_tools.plot_hist(data1, bins=50, log_scale=True)
#         plot_tools.set_labels("X-axis", "Frequency")
#         plot_tools.set_limits(xlim=(-3, 3))

#     plot_first_hist(plot_tools)
#     plot_tools.next_plot()

#     @plot_tools.with_title("2D Histogram")
#     @plot_tools.with_grid(color='blue')
#     def plot_second_hist(plot_tools):
#         x2 = np.random.normal(size=1000)
#         y2 = np.random.normal(size=1000)
#         plot_tools.plot_hist2d(x2, y2, bins=50, log_scale=True)
#         plot_tools.set_labels("X2-axis", "Y2-axis")
#         plot_tools.set_limits(xlim=(-3, 3), ylim=(-3, 3))

#     plot_second_hist(plot_tools)

#     # Save the figure to a PDF file
#     plot_tools.save_to_pdf("plots.pdf")
