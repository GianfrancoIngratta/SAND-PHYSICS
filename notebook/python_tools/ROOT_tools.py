import ROOT
import time

class ROOT_tools:

    def __init__(self, arg_file = ""):
        self.file = arg_file
    
    def _check_TH1D_(self, h: ROOT.TH1D):
        if not isinstance(h, ROOT.TH1D) or h.GetEntries() == 0:
            raise ValueError("input hist has to be a non empty ROOT.TH1D")
        else:
            return True

    def FillTH1D(self, iterable: list, histogram_name: str, title: str, nbins: int , x_min: int, x_max: int) -> ROOT.TH1D:
        # fill TH1D from python-like vector
        # prevent memory leak
        if len(iterable)==0: raise ValueError("input vector is empty")
        histogram_name = histogram_name + " time:" + str(int(time.time()))
        histogram = ROOT.TH1D(histogram_name, title, nbins, x_min, x_max)
        for entry in iterable: histogram.Fill(entry)
        return histogram

    def FitTH1D_w_gauss(self, hist: ROOT.TH1D, gauss_range: tuple, fit_range: tuple) -> tuple:
        # fit TH1D with gauss function and return hist, mean and sigma of the fit
        self._check_TH1D_(hist)
        gaussian_func = ROOT.TF1("gaussian_func", "gaus", gauss_range[0], gauss_range[1])
        gaussian_func.SetRange(fit_range[0], fit_range[1])
        hist.Fit(gaussian_func, "R")
        mean = gaussian_func.GetParameter(1)
        sigma = gaussian_func.GetParameter(2)
        return (hist, mean, sigma)
    
    def FitTH1D_w_chi2(self, hist: ROOT.TH1D, chi2_range: tuple, fit_range: tuple) -> tuple:
        # fit TH1D with chi2 function and return hist, mean and sigma of the fit
        self._check_TH1D_(hist)
        chi2_func = ROOT.TF1("chi2_func", "[0]*ROOT::Math::chisquared_pdf(x, [1])", chi2_range[0], chi2_range[1])
        # Initial guess for the parameters: scale and degrees of freedom
        chi2_func.SetParameters(1, 1)
        # Set the range for the fit
        chi2_func.SetRange(fit_range[0], fit_range[1])
        # Fit the histogram with the chi-squared function
        hist.Fit(chi2_func, "R")
        # Extract the parameters: scale and degrees of freedom
        scale = chi2_func.GetParameter(0)
        dof = chi2_func.GetParameter(1)
        return (hist, scale, dof)
    
    def PlotTH1D(self, hist: ROOT.TH1D, canvas_name: str, canvas_dimensions = (500,800)) -> None:
        self._check_TH1D_(hist)
        canvas_name = "{} time: {}".format(canvas_name, str(int(time.time()))[-4:])
        ROOT.gStyle.SetOptFit(1011)
        canvas = ROOT.TCanvas(canvas_name, "Canvas", canvas_dimensions[0], canvas_dimensions[1])
        hist.Draw()
        canvas.Draw()
    