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

default_pt_categories = [[450., 500.], 	[500., 550.], [550., 600.], [600., 675.], [675., 800.], [800., 1000.]]

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

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description="Make simple datacards/workspaces/run scripts/limit plots")
	parser.add_argument('--make_datacards', action="store_true", help="Make datacards and workspaces")
	parser.add_argument('--make_plots', action="store_true", help="Make limit plots")
	parser.add_argument('--nrho', type=int, default=2, help="Degree of rho polynomial")
	parser.add_argument('--npt', type=int, default=1, help="Degree of pt polynomial")
	args = parser.parse_args()

	if args.make_datacards:
		histogram_files = {
			"vbf":"/uscms/home/kkwok/work/Hbb/CMSSW_8_1_0/src/ZPrimePlusJet/fitting/PbbJet/QG/QGquark/hist_1DZbb_pt_scalesmear.root",
			"ggf":"/uscms/home/kkwok/work/Hbb/CMSSW_8_1_0/src/ZPrimePlusJet/fitting/PbbJet/QG/QGgluon/hist_1DZbb_pt_scalesmear.root"
		}
		datacard_dir = "/uscms/home/dryu/HbbPlusJet/data/datacards/"
		os.system("mkdir -pv {}".format(datacard_dir))
		project = CombineProject(datacard_dir)
		xvar = project.create_xvar("msd", [40., 201.])

		# Get raw histograms
		h2d = {}
		h2d_slices = {}
		for tag in ["ggf", "vbf"]:
			histogram_file = TFile(histogram_files[tag], "READ")
			for box in ["pass", "fail"]:
				region = "{}_{}".format(tag, box)
				h2d[region] = {}
				h2d_slices[region] = {}

				# Load 2D histograms (ex.: ZPrime150_failn2_faildbtag_pt_vs_msd)
				for process in backgrounds + signals:
					hname = "{}_{}".format(process, box)
					if process in signals + ["zqq", "wqq"]:
						hname += "_matched"
					h2d[region][process] = histogram_file.Get(hname)
					if not h2d[region][process]:
						print "[rhalphabet_limits] ERROR : Couldn't load histogram {} from file {}".format("{}_{}_{}".format(process, tag, "pt_vs_msd"), histogram_file.GetPath())
						sys.exit(1)
					h2d[region][process].SetName("{}_{}".format(process, region))
					h2d[region][process].SetDirectory(0)

				# Normalize QCD to match data normalization
				#h2d_data_tmp = histogram_file.Get("data_obs_{}_pt_vs_msd".format(tag))
				#int_data = h2d_data_tmp.Integral()
				#int_bkgds = {}
				#for process in backgrounds:
				#	int_bkgds[process] = h2d[tag][process].Integral()
				#if h2d[tag]["qcd"].Integral() == 0:
				#	print "[rhalphabet_limits] ERROR : QCD histogram has zero integral. File={}, hist={}".format(histogram_file.GetPath(), "{}_{}_{}".format("qcd", tag, "pt_vs_msd"))
				#h2d[tag]["qcd"].Scale((int_data - sum([int_bkgds[process] for process in backgrounds if process != "qcd"])) / h2d[tag]["qcd"].Integral())

				# Slice 2D histograms
				for process in backgrounds + signals:
					h2d_slices[region][process] = SliceSignalHistogram(h2d[region][process], pt_categories=default_pt_categories)

				# Data histogram = pseudodata from SM
				for process in backgrounds + signals:
					if not "data_obs" in h2d[region]:
						h2d[region]["data_obs"] = h2d[region][process].Clone()
						h2d[region]["data_obs"].SetName("data_obs_{}".format(region))
						h2d[region]["data_obs"].SetDirectory(0)
					else:
						h2d[region]["data_obs"].Add(h2d[region][process])
				h2d_slices[region]["data_obs"] = SliceSignalHistogram(h2d[region]["data_obs"], pt_categories=default_pt_categories)
			# End loop over tags
		
		# Apply rho boundaries, and set data to zero if the prediction is zero
		rho_range = [-6.0, -2.1]

		for pt_index in xrange(len(default_pt_categories)):
			for tag in ["ggf", "vbf"]:
				pass_tag = "{}_pass".format(tag)
				fail_tag = "{}_fail".format(tag)

				for pf_tag in [pass_tag, fail_tag]:
					for xbin in xrange(1, h2d_slices[pf_tag]["data_obs"][pt_index].GetNbinsX()+1):
						bin_mass = h2d_slices[pf_tag]["data_obs"][pt_index].GetXaxis().GetBinCenter(xbin)
						bin_pt = default_pt_categories[pt_index][0] + (default_pt_categories[pt_index][1]-default_pt_categories[pt_index][0]) * 0.3
						bin_rho = 2. * math.log(bin_mass / bin_pt)
						if bin_rho < rho_range[0] or bin_rho > rho_range[1]:
							for process in backgrounds + signals + ["data_obs"]:
								h2d_slices[pf_tag][process][pt_index].SetBinContent(xbin, 0.)
								h2d_slices[pf_tag][process][pt_index].SetBinError(xbin, 0.)
				# End loop over N2 pass/fail

				# If the fail bin 0 entries, set pass=0.
				if h2d_slices[fail_tag]["data_obs"][pt_index].GetBinContent(xbin) == 0:
					h2d_slices[fail_tag]["data_obs"][pt_index].SetBinContent(xbin, 1.e-10)
					h2d_slices[pass_tag]["data_obs"][pt_index].SetBinContent(xbin, 0)
					h2d_slices[pass_tag]["data_obs"][pt_index].SetBinError(xbin, 0)
			# End loop over ggf/vbf
		# End loop over pt categories

		
		# Reorganize maps: <category name, <process, histogram>>, where category namee = e.g. ptcatN, or dbtagpass_ptcatN
		region_names = []
		pass_histograms = {}
		fail_histograms = {}
		region_pts = {}
		region_cats = {}
		for ptcat in xrange(1, len(default_pt_categories) + 1):
			for tag in ["ggf", "vbf"]:
				region_name = "{}_ptcat{}".format(tag, ptcat)
				region_names.append(region_name)
				region_pts[region_name] = default_pt_categories[ptcat-1]
				region_cats[region_name] = ptcat
				pass_histograms[region_name] = {}
				fail_histograms[region_name] = {}
				for process in backgrounds + signals + ["data_obs"]:
					pass_histograms[region_name][process] = h2d_slices["{}_pass".format(tag)][process][ptcat-1]
					fail_histograms[region_name][process] = h2d_slices["{}_fail".format(tag)][process][ptcat-1]

		#CreateCombineProject(region_names, pass_histograms, fail_histograms, jet_type=jet_type, datacard_directory=datacard_dir, signal_name=signal, data_name="data_obs", background_names=backgrounds, region_pts=region_pts)
		project = CombineProject(datacard_dir)	
		# Construct signal regions
		region_containers = {}
		for region_name in region_names:
			region_containers[region_name] = RhalphabetRegion(region_name)
			region_containers[region_name].add_xvar(xvar)

			# def add_data(self, data_name, pass_hist, fail_hist):
			region_containers[region_name].add_data("data_obs", pass_histograms[region_name]["data_obs"], fail_histograms[region_name]["data_obs"])
			# def add_signal(self, signal_name, pass_hist, fail_hist, normalization_var=None):

			for signal in signals:
				region_containers[region_name].add_signal(signal, pass_histograms[region_name][signal], fail_histograms[region_name][signal])

			for background_name in backgrounds:
				if background_name == "qcd":
					continue
				else:
					# def add_simple_background(self, bkgd_name, pass_hist, fail_hist, normalization_var=None):
					region_containers[region_name].add_simple_background(background_name, pass_histograms[region_name][background_name], fail_histograms[region_name][background_name])
			# y_value = pT 1/3 of the way up the bin
			y_value = (1. * region_pts[region_name][0] + 2. * region_pts[region_name][1]) / 3.
			rh_tag = ""
			if "ggf" in region_name:
				rh_tag = "_ggf"
			elif "vbf" in region_name:
				rh_tag = "_vbf"
			region_containers[region_name].add_rhbackground("qcd", y_value, n_x=args.nrho, n_y=args.npt, x_range=[-7., 0.], y_range=[400., 1000.], rh_tag=rh_tag)

			###################
			### Systematics ###
			###################
			# def add_norm_systematic(self, syst_name, process, unc_value):
			#systematics_norm = systematics.GetNormalizationSystematics(args.jet_type, args.wp)
			#if "pass_dbtag" in region_name:
			#	pass_tag = "passn2_passdbtag"
			#	fail_tag = "failn2_passdbtag"
			#else:
			#	pass_tag = "passn2_faildbtag"
			#	fail_tag = "failn2_faildbtag"
			#for process in backgrounds + [signal]:
			#	# Skip QCD, which is handled independently
			#	if process == "qcd":
			#		continue
			#	if process in systematics_norm["lumi"]:
			#		region_containers[region_name].add_norm_systematic("lumi", process, (systematics_norm["lumi"][process], systematics_norm["lumi"][process]))
			#	if process in systematics_norm["muveto"]:
			#		region_containers[region_name].add_norm_systematic("muveto", process, (systematics_norm["muveto"][process], systematics_norm["muveto"][process]))
			#	if process in systematics_norm["eleveto"]:
			#		region_containers[region_name].add_norm_systematic("eleveto", process, (systematics_norm["eleveto"][process], systematics_norm["eleveto"][process]))

			#	if process in systematics_norm["veff"]:
			#		region_containers[region_name].add_norm_systematic("veff", process, (systematics_norm["veff"][process][pass_tag], systematics_norm["veff"][process][fail_tag]))
			#	if process in systematics_norm["bbeff"]:
			#		region_containers[region_name].add_norm_systematic("bbeff", process, (systematics_norm["bbeff"][process][pass_tag], systematics_norm["bbeff"][process][fail_tag]))

			#	if process in systematics_norm["znormQ"]:
			#		region_containers[region_name].add_norm_systematic("znormQ", process, (systematics_norm["znormQ"][process], systematics_norm["znormQ"][process]))
			#	if process in systematics_norm["znormEW"]:
			#		region_containers[region_name].add_norm_systematic("znormEW", process, (systematics_norm["znormEW"][process][region_cats[region_name]-1], systematics_norm["znormEW"][process][region_cats[region_name]-1]))
			#	if process in systematics_norm["wznormEW"]:
			#		region_containers[region_name].add_norm_systematic("wznormEW", process, (systematics_norm["wznormEW"][process][region_cats[region_name]-1], systematics_norm["wznormEW"][process][region_cats[region_name]-1]))

			#	if process in systematics_norm["hpt"]:
			#		region_containers[region_name].add_norm_systematic("hpt", process, (systematics_norm["hpt"][process], systematics_norm["hpt"][process]))

			#	if process in systematics_norm["JER"]:
			#		region_containers[region_name].add_norm_systematic("JER", process, (systematics_norm["JER"][process][pass_tag], systematics_norm["JER"][process][fail_tag]))
			#	if process in systematics_norm["JES"]:
			#		region_containers[region_name].add_norm_systematic("JES", process, (systematics_norm["JES"][process][pass_tag], systematics_norm["JES"][process][fail_tag]))
			#	if process in systematics_norm["PU"]:
			#		region_containers[region_name].add_norm_systematic("PU", process, (systematics_norm["PU"][process][pass_tag], systematics_norm["PU"][process][fail_tag]))
			#	if process in systematics_norm["Trigger"]:
			#		region_containers[region_name].add_norm_systematic("Trigger", process, (systematics_norm["Trigger"][process][pass_tag], systematics_norm["Trigger"][process][fail_tag]))

			project.add_region(region_containers[region_name])
		project.write()

		# Run all script
		with open("{}/run_limits.sh".format(datacard_dir), 'w') as run_all_script:
			run_all_script.write("#!/bin/bash\n")
			run_all_script.write("cd {}\n".format(datacard_dir))
			run_all_script.write("combine -M AsymptoticLimits datacard_total.txt >& limits.log\n")
			run_all_script.write("cd -\n".format(datacard_dir))
			run_all_script.close()


	if args.make_plots:
		from DAZSLE.PhiBBPlusJet.limit_plot import LimitPlot
		import array
		signal_model_masses = {
			"Sbb":array.array('d', [50,100,125,200,300,350,400,500]),
			"ZPrime":array.array('d', [75, 125, 150, 175, 225, 250, 300, 400])
		}

		for signal_model in ["Sbb", "ZPrime"]:
			limit_whats = ["-2exp", "-1exp", "exp", "+1exp", "+2exp", "obs"]
			limits_btag = {}
			limits_nobtag = {}
			for limit_what in limit_whats:
				limits_btag[limit_what] = array.array('d', [])
				limits_nobtag[limit_what] = array.array('d', [])
			print "{} - {}".format(args.jet_type, signal_model)
			for mass in signal_model_masses[signal_model]:
				limit_file_btag = TFile("/uscms/home/dryu/HbbPlusJet/data/datacards/simple/{}/{}{}/higgsCombineTest.AsymptoticLimits.mH120.root".format(args.jet_type, signal_model, int(mass)), "READ")
				print "Opening {}".format(limit_file_btag.GetPath())
				t_btag = limit_file_btag.Get("limit")
				for i, limit_what in enumerate(limit_whats):
					t_btag.GetEntry(i)
					limits_btag[limit_what].append(t_btag.GetLeaf("limit").GetValue(0))

				limit_file_nobtag = TFile("/uscms/home/dryu/HbbPlusJet/data/datacards/simple_nobtag/{}/{}{}/higgsCombineTest.AsymptoticLimits.mH120.root".format(args.jet_type, signal_model, int(mass)), "READ")
				t_nobtag = limit_file_nobtag.Get("limit")
				for i, limit_what in enumerate(limit_whats):
					t_nobtag.GetEntry(i)
					limits_nobtag[limit_what].append(t_nobtag.GetLeaf("limit").GetValue(0))
			# End loop over masses
			 
			print "Exp limits with b-tag: ",
			print limits_btag["exp"]
			print "Exp limits without b-tag: ",
			print limits_nobtag["exp"]	
			
			# Make TGraphs
			limit_tgraphs_btag = {}
			limit_tgraphs_nobtag = {}
			for limit_what in limit_whats:
				limit_tgraphs_btag[limit_what] = TGraph(len(signal_model_masses[signal_model]), signal_model_masses[signal_model], limits_btag[limit_what])
				limit_tgraphs_nobtag[limit_what] = TGraph(len(signal_model_masses[signal_model]), signal_model_masses[signal_model], limits_nobtag[limit_what])

			# Use the LimitPlot class to make the fill objects
			limit_plot_btag = LimitPlot()
			limit_plot_btag.load_limit_graph(limit_tgraphs_btag["obs"])
			limit_plot_btag.load_limit_graph_exp(-2, limit_tgraphs_btag["-2exp"])
			limit_plot_btag.load_limit_graph_exp(-1, limit_tgraphs_btag["-1exp"])
			limit_plot_btag.load_limit_graph_exp(0, limit_tgraphs_btag["exp"])
			limit_plot_btag.load_limit_graph_exp(1, limit_tgraphs_btag["+1exp"])
			limit_plot_btag.load_limit_graph_exp(2, limit_tgraphs_btag["+2exp"])
			limit_plot_btag.draw("limit_plot_btag")
			limit_fills_btag_2sig = limit_plot_btag._limit_fills_exp[2]
			limit_fills_btag_1sig = limit_plot_btag._limit_fills_exp[1]

			limit_plot_nobtag = LimitPlot()
			limit_plot_nobtag.load_limit_graph(limit_tgraphs_nobtag["obs"])
			limit_plot_nobtag.load_limit_graph_exp(-2, limit_tgraphs_nobtag["-2exp"])
			limit_plot_nobtag.load_limit_graph_exp(-1, limit_tgraphs_nobtag["-1exp"])
			limit_plot_nobtag.load_limit_graph_exp(0, limit_tgraphs_nobtag["exp"])
			limit_plot_nobtag.load_limit_graph_exp(1, limit_tgraphs_nobtag["+1exp"])
			limit_plot_nobtag.load_limit_graph_exp(2, limit_tgraphs_nobtag["+2exp"])
			limit_plot_nobtag.draw("limit_plot_nobtag")
			limit_fills_nobtag_2sig = limit_plot_nobtag._limit_fills_exp[2]
			limit_fills_nobtag_1sig = limit_plot_nobtag._limit_fills_exp[1]

			# Drawing
			yellow_transparent = TColor.GetColorTransparent(kOrange, 0.5)
			green_transparent = TColor.GetColorTransparent(kGreen, 0.5)
			limit_fills_btag_2sig.SetFillColor(yellow_transparent)
			limit_fills_btag_1sig.SetFillColor(green_transparent)
			limit_fills_nobtag_2sig.SetFillColor(yellow_transparent)
			limit_fills_nobtag_1sig.SetFillColor(green_transparent)

			c = TCanvas("c_limitcomparison_{}_{}".format(signal_model, args.jet_type), "Limit comparison", 800, 600)
			frame = TH1D("frame", "frame", 100, 0., 500.)
			if signal_model == "Sbb":
				frame.SetMinimum(0.02)
				frame.SetMaximum(10.)
			else:
				frame.SetMinimum(1.)
				frame.SetMaximum(100.)
			frame.GetXaxis().SetTitle("m_{X} [GeV]")
			frame.GetYaxis().SetTitle("#sigma#timesBR(jj)#times A #times #epsilon [pb]")
			c.SetLogy()
			frame.Draw()
			limit_fills_btag_2sig.Draw("f")
			limit_fills_btag_1sig.Draw("f")
			limit_fills_nobtag_2sig.Draw("f")
			limit_fills_nobtag_1sig.Draw("f")

			limit_tgraphs_btag["exp"].SetLineColor(kBlack)
			limit_tgraphs_btag["exp"].SetLineStyle(3)
			limit_tgraphs_btag["exp"].SetLineWidth(2)
			limit_tgraphs_btag["exp"].SetMarkerStyle(20)
			limit_tgraphs_btag["exp"].SetMarkerSize(0)
			limit_tgraphs_btag["exp"].Draw("l")

			limit_tgraphs_nobtag["exp"].SetLineColor(kBlue)
			limit_tgraphs_nobtag["exp"].SetLineStyle(3)
			limit_tgraphs_nobtag["exp"].SetLineWidth(2)
			limit_tgraphs_nobtag["exp"].SetMarkerStyle(20)
			limit_tgraphs_nobtag["exp"].SetMarkerSize(0)
			limit_tgraphs_nobtag["exp"].Draw("l")

			limit_tgraphs_btag["obs"].SetLineColor(kBlack)
			limit_tgraphs_btag["obs"].SetLineStyle(1)
			limit_tgraphs_btag["obs"].SetLineWidth(2)
			limit_tgraphs_btag["obs"].SetMarkerStyle(20)
			limit_tgraphs_btag["obs"].SetMarkerSize(0)
			#limit_tgraphs_btag["obs"].Draw("l")

			limit_tgraphs_nobtag["obs"].SetLineColor(kBlue)
			limit_tgraphs_nobtag["obs"].SetLineStyle(1)
			limit_tgraphs_nobtag["obs"].SetLineWidth(2)
			limit_tgraphs_nobtag["obs"].SetMarkerStyle(20)
			limit_tgraphs_nobtag["obs"].SetMarkerSize(0)
			#limit_tgraphs_nobtag["obs"].Draw("l")

			l = TLegend(0.23, 0.6, 0.4, 0.8)
			l.SetFillStyle(0)
			l.SetBorderSize(0)
			l.SetHeader("95% CL upper limits")
			l.AddEntry(limit_tgraphs_btag["exp"], "Exp. limit (w/b-tag)", "l")
			l.AddEntry(limit_tgraphs_nobtag["exp"], "Exp. limit (wo/b-tag)", "l")
			#l.AddEntry(limit_tgraphs_btag["exp"], "Obs. limit (w/b-tag)", "l")
			#l.AddEntry(limit_tgraphs_nobtag["exp"], "Obs. limit (wo/b-tag)", "l")
			l.Draw()

			c.SaveAs("/uscms/home/dryu/HbbPlusJet/data/limits/figures/{}.pdf".format(c.GetName()))
		# End loop over signal models
