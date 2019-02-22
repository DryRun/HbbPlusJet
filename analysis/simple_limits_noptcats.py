import os
import sys
from DAZSLE.DAZSLECommon.combine_helper.combine_project import CombineProject
from DAZSLE.DAZSLECommon.combine_helper.region import Region
from DAZSLE.DAZSLECommon.combine_helper.rhalphabet_region import RhalphabetRegion
from ROOT import TFile, TGraph, TColor, gROOT, kOrange, kGreen, TCanvas, TLegend, TH1D, kBlack, kBlue
gROOT.SetBatch(True)
import math

#import systematics

# Load histograms
backgrounds = ["qcd", "tqq", "zqq", "wqq"]
signals = ["hqq125", "tthqq125", "vbfhqq125", "whqq125", "zhqq125"]

default_pt_categories = [[450., 1000.]]

# Slice up the pt vs msd histogram into 1D msd histograms
def SliceSignalHistogram(hist2d, pt_categories=default_pt_categories):
	# Make sure the requested pt bins correspond to boundaries
	pt_axis = hist2d.GetYaxis()
	histogram_pt_boundaries = []
	for bin in xrange(1, pt_axis.GetNbins() + 1):
		histogram_pt_boundaries.append(pt_axis.GetBinLowEdge(bin))
	histogram_pt_boundaries.append(pt_axis.GetBinUpEdge(pt_axis.GetNbins()))

	for pt_category in pt_categories:
		if not pt_category[0] in histogram_pt_boundaries:
			print "[run_histograms::SliceSignalHistogram] ERROR : Bin boundary {} does not correspond to a histogram bin boundary.".format(pt_category[0])
			print histogram_pt_boundaries
			sys.exit(1)
		if not pt_category[1] in histogram_pt_boundaries:
			print "[run_histograms::SliceSignalHistogram] ERROR : Bin boundary {} does not correspond to a histogram bin boundary.".format(pt_category[1])
			print histogram_pt_boundaries
			sys.exit(1)

	histogram_slices = {}
	for islice in xrange(len(pt_categories)):
		ptmin = pt_categories[islice][0]
		ptmax = pt_categories[islice][1]
		binmin = 1e10
		binmax = -1e10
		for bin in xrange(1, hist2d.GetNbinsY() + 1):
			low_edge = hist2d.GetYaxis().GetBinLowEdge(bin)
			up_edge = hist2d.GetYaxis().GetBinUpEdge(bin)
			# Is this bin inside this pt slice (+epsilon)?
			if ptmin - 1.e-5 < low_edge and up_edge < ptmax + 1.e-5:
				if bin < binmin:
					binmin = bin
				if bin > binmax:
					binmax = bin
		histogram_slices[islice] = hist2d.ProjectionX(hist2d.GetName() + "_" + str(islice), binmin, binmax)
		histogram_slices[islice].SetDirectory(0)
	return histogram_slices

def MassScaling(hist_nominal):
	# Measured mean of W in data and MC
	mean_data = 82.657
	mean_data_err = 0.313
	mean_mc = 82.548
	mean_mc_err = 0.191

	# Scale factors to transform MC shape into data shape
	mass_shift = mean_data / mean_mc
	mass_shift_unc = math.sqrt((mean_data_err / mean_data) * (mean_data_err / mean_data) + (mean_mc_err / mean_mc) * (mean_mc_err / mean_mc)) * 10.  # (10 sigma shift)

	# Mass shift histogram
	shift_param = mass * (1. -mass_shift)
	xvar =   RooRealVar("x", "x", hist_nominal.GetXaxis().GetXmin(), hist_nominal.GetXaxis().GetXmax())
	xvar.setBins(hist_nominal.GetNbinsX())
	dxvar = RooRealVar("dx", "dx", 0., -10., 10.)
	shiftvar = RooFormulaVar("shift", "x-dx", RooArgList(xvar, dxvar))

	datahist_shift = RooDataHist(hist_nominal.GetName()+"_rdh_shift" ,hist_nominal.GetName()+"_rdh_shift" ,r.RooArgList(xvar), hist_nominal)
	histpdf_shift = RooHistPdf(hist_nominal.GetName() + "_rhp_shift", hist_nominal.GetName() + "_rhp_shift", RooArgList(shiftvar), RooArgList(xvar), datahist_shift)

	shiftvar.setVal(shift_param)
	hist_massscale_up = histpdf_shift.createHistogram("x")
	hist_massscale_up.SetTitle(hist_nominal.GetTitle() + "_scaleUp")
	hist_massscale_up.SetName(hist_nominal.GetName() + "_scaleUp")
	hist_massscale_up.Scale(hist_nominal.Integral())

	shiftvar.setVal(-1. * shift_param)
	hist_massscale_down = histpdf_shift.createHistogram("x")
	hist_massscale_down.SetTitle(hist_nominal.GetTitle() + "_scaleDown")
	hist_massscale_down.SetName(hist_nominal.GetName() + "_scaleDown")
	hist_massscale_down.Scale(hist_nominal.Integral())

	return (hist_massscale_up, hist_massscale_down)


def MassSmearing(hist_nominal):
	# Measured sigma of W in data and MC
	sigma_data = 8.701
	sigma_data_err = 0.433
	sigma_mc = 8.027
	sigma_mc_err = 0.607

	# Scale factors to transform MC shape into data shape
	res_shift = s_data / s_mc
	res_shift_unc = math.sqrt((sigma_data_err / sigma_data) * (sigma_data_err / sigma_data) + (sigma_mc_err / sigma_mc) * (sigma_mc_err / sigma_mc)) * 2.  # (2 sigma shift)

	# Mass smear histogram
	smear_param = res_shift - 1
	xvar =   RooRealVar("x", "x", hist_nominal.GetXaxis().GetXmin(), hist_nominal.GetXaxis().GetXmax())
	xvar.setBins(hist_nominal.GetNbinsX())
	dxvar = RooRealVar("dx", "dx", 0., 0., 2.)
	smearvar = RooFormulaVar("smear","(x-{})/dx+{}".format(hist_nominal.GetMean(), hist_nominal.GetMean()), RooArgList(xvar, dxvar))

	datahist_smear = RooDataHist(hist_nominal.GetName()+"_rdh_smear" ,hist_nominal.GetName()+"_rdh_smear" ,r.RooArgList(xvar), hist_nominal)
	histpdf_smear = RooHistPdf(hist_nominal.GetName() + "_rhp_smear", hist_nominal.GetName() + "_rhp_smear", RooArgList(smearvar), RooArgList(xvar), datahist_smear)

	smearvar.setVal(1. + smear_param)
	hist_massres_up = histpdf_smear.createHistogram("x")
	hist_massres_up.SetTitle(hist_nominal.GetTitle() + "_smearUp")
	hist_massres_up.SetName(hist_nominal.GetName() + "_smearUp")
	hist_massres_up.Scale(hist_nominal.Integral())

	smearvar.setVal(1. - smear_param)
	hist_massres_down = histpdf_smear.createHistogram("x")
	hist_massres_down.SetTitle(hist_nominal.GetTitle() + "_smearDown")
	hist_massres_down.SetName(hist_nominal.GetName() + "_smearDown")
	hist_massres_down.Scale(hist_nominal.Integral())

	return (hist_massres_up, hist_massres_down)



if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Make simple datacards/workspaces/run scripts/limit plots")
	parser.add_argument('--make_datacards', action="store_true", help="Make datacards and workspaces")
	parser.add_argument('--make_plots', action="store_true", help="Make limit plots")
	args = parser.parse_args()

	if args.make_datacards:
		histogram_files = {
			"vbf":"/uscms/home/kkwok/work/Hbb/CMSSW_8_1_0/src/ZPrimePlusJet/fitting/PbbJet/QG/QGquark/hist_1DZbb_pt_scalesmear.root",
			"ggf":"/uscms/home/kkwok/work/Hbb/CMSSW_8_1_0/src/ZPrimePlusJet/fitting/PbbJet/QG/QGgluon/hist_1DZbb_pt_scalesmear.root"
		}
		datacard_dir = "/uscms/home/dryu/HbbPlusJet/data/datacards/simple_noptcats"
		os.system("mkdir -pv {}".format(datacard_dir))
		project = CombineProject(datacard_dir)
		xvar = project.create_xvar("msd", [40., 201.])

		# Get raw histograms
		h2d = {}
		h2d_slices = {}
		h2d_syst = {}
		h2d_slices_syst = {}
		for tag in ["ggf", "vbf"]:
			histogram_file = TFile(histogram_files[tag], "READ")
			for box in ["pass"]:
				region = "{}_{}".format(tag, box)

				# Load 2D histograms (ex.: ZPrime150_failn2_faildbtag_pt_vs_msd)
				for process in backgrounds + signals:
					hname = "{}_{}".format(process, box)
					if process in signals + ["zqq", "wqq"]:
						hname += "_matched"
					if not process in h2d:
						h2d[process] = histogram_file.Get(hname)
						h2d[process].SetName("{}".format(process))
						h2d[process].SetDirectory(0)
					else:
						h2d[process].Add(histogram_file.Get(hname))

					h2d_syst[process] = {}
					h2d_slices_syst[process] = {}
					for syst in ["JER", "JES", "Pu", "trigger"]:
						for direction in ["Up", "Down"]:
							hname = "{}_{}_{}{}".format(process, box, syst, direction)
							if not syst + direction in h2d_syst[process]:
								h2d_syst[process][syst + direction] = histogram_file.Get(hname)
								h2d_syst[process][syst + direction].SetName("{}_{}{}".format(process, syst, direction))
								h2d_syst[process][syst + direction].SetDirectory(0)
							else:
								h2d_syst[process][syst + direction].Add(histogram_file.Get(hname))

				# Normalize QCD to match data normalization
				#h2d_data_tmp = histogram_file.Get("data_obs_{}_pt_vs_msd".format(tag))
				#int_data = h2d_data_tmp.Integral()
				#int_bkgds = {}
				#for process in backgrounds:
				#	int_bkgds[process] = h2d[tag][process].Integral()
				#if h2d[tag]["qcd"].Integral() == 0:
				#	print "[rhalphabet_limits] ERROR : QCD histogram has zero integral. File={}, hist={}".format(histogram_file.GetPath(), "{}_{}_{}".format("qcd", tag, "pt_vs_msd"))
				#h2d[tag]["qcd"].Scale((int_data - sum([int_bkgds[process] for process in backgrounds if process != "qcd"])) / h2d[tag]["qcd"].Integral())
			# End loop over pass/fail
		# End loop over ggf/vbf
		
		# Slice 2D histograms
		for process in backgrounds + signals:
			h2d_slices[process] = SliceSignalHistogram(h2d[process], pt_categories=default_pt_categories)
			for syst in ["JER", "JES", "Pu", "trigger"]:
				for direction in ["Up", "Down"]:
					h2d_slices_syst[process][syst + direction] = SliceSignalHistogram(h2d_syst[process][syst + direction], pt_categories=default_pt_categories)

		# Data histogram = pseudodata from SM
		for process in backgrounds + signals:
			if not "data_obs" in h2d:
				h2d["data_obs"] = h2d[process].Clone()
				h2d["data_obs"].SetName("data_obs")
				h2d["data_obs"].SetDirectory(0)
			else:
				h2d["data_obs"].Add(h2d[process])
		h2d_slices["data_obs"] = SliceSignalHistogram(h2d["data_obs"], pt_categories=default_pt_categories)
		

		# Apply rho boundaries, and set data to zero if the prediction is zero
		rho_range = [-6.0, -2.1]

		for pt_index in xrange(len(default_pt_categories)):
			for xbin in xrange(1, h2d_slices["data_obs"][pt_index].GetNbinsX()+1):
				bin_mass = h2d_slices["data_obs"][pt_index].GetXaxis().GetBinCenter(xbin)
				bin_pt = default_pt_categories[pt_index][0] + (default_pt_categories[pt_index][1]-default_pt_categories[pt_index][0]) * 0.3
				bin_rho = 2. * math.log(bin_mass / bin_pt)
				if bin_rho < rho_range[0] or bin_rho > rho_range[1]:
					for process in backgrounds + signals + ["data_obs"]:
						h2d_slices[process][pt_index].SetBinContent(xbin, 0.)
						h2d_slices[process][pt_index].SetBinError(xbin, 0.)
		
		#CreateCombineProject(region_names, pass_histograms, fail_histograms, jet_type=jet_type, datacard_directory=datacard_dir, signal_name=signal, data_name="data_obs", background_names=backgrounds, region_pts=region_pts)
		project = CombineProject(datacard_dir)	
		# Construct signal region
		signal_region = Region("ptcat0")
		signal_region.add_xvar(xvar)

		# def add_data(self, data_name, pass_hist, fail_hist):
		signal_region.add_data("data_obs", h2d_slices["data_obs"][0])
		# def add_signal(self, signal_name, pass_hist, fail_hist, normalization_var=None):

		for signal in signals:
			signal_region.add_signal(signal, h2d_slices[signal][0])

		for background_name in backgrounds:
			signal_region.add_background(background_name, h2d_slices[background_name][0])

		###################
		### Systematics ###
		###################
		# def add_norm_systematic(self, syst_name, process, unc_value=None, unc_pass_fail=None):
		
		# Luminosity
		for process in backgrounds + signals:
			signal_region.add_norm_systematic("lumi", process, unc_value=1.025)

		# Higgs pT
		signal_region.add_norm_systematic("hqq125pt", "hqq125", unc_value=1.3)

		# N2 and double b-tag SFs
		for process in backgrounds + signals:
			if process in ["qcd", "tqq"]:
				continue
			signal_region.add_norm_systematic("veff", process, unc_value=1.043)

			# beff: need to derive fail region SF from pass+fail normalization constraint
			bbsf = 0.91
			dbbsf = 0.03
			signal_region.add_norm_systematic("beff", process, unc_value=1 + dbbsf)

		# Electroweak normalizations
		for process in ["wqq", "zqq"]:
			signal_region.add_norm_systematic("znormQ", process, unc_value=1.1)
			signal_region.add_norm_systematic("wznormEW", process, unc_value=1.05)
			signal_region.add_norm_systematic("znormEW", process, unc_value=1.15)

		# JES, JER, Pu, trigger
		for systematic in ["JER", "JES", "Pu", "trigger"]:
			for process in backgrounds + signals:
				if process == "qcd":
					continue

				norm = h2d_slices[background_name][0].Integral()
				norm_up = h2d_slices_syst[background_name][syst + "Up"][0].Integral()
				norm_down = h2d_slices_syst[background_name][syst + "Down"][0].Integral()
				syst_value = 1.0 + (abs(norm_up - norm) + abs(norm - norm_down)) / (2. * norm)
				signal_region.add_norm_systematic(systematic, process, unc_value=syst_value)

		# Lepton veto
		for process in backgrounds + signals:
			if process == "qcd":
				continue
			signal_region.add_norm_systematic("muveto", process, unc_value=1.005)
			signal_region.add_norm_systematic("eleveto", process, unc_value=1.005)

		# Shape systematics 
		# Mass scale and resolution
		#for process in ["wqq", "zqq"] + signals:
		#	hists_massscale_pass = MassScaling(pass_histograms[region_name][process])
		#	hists_massscale_fail = MassScaling(fail_histograms[region_name][process])
		#	signal_region.add_shape_systematic("massscale", process, hists_massscale_pass, hists_massscale_fail)

		#	hists_masssmear_pass = MassSmearing(pass_histograms[region_name][process])
		#	hists_masssmear_fail = MassSmearing(fail_histograms[region_name][process])
		#	signal_region.add_shape_systematic("masssmear", process, hists_masssmear_pass, hists_masssmear_fail)


		project.add_region(signal_region)
		project.write()

		# Run all script
		with open("{}/run_limits.sh".format(datacard_dir), 'w') as run_all_script:
			run_all_script.write("#!/bin/bash\n")
			run_all_script.write("cd {}\n".format(datacard_dir))
			run_all_script.write("combine -M AsymptoticLimits datacard_total.txt >& limits.log\n")
			run_all_script.write("cd -\n".format(datacard_dir))
			run_all_script.close()


