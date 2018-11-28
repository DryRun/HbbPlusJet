import os
import sys
from DAZSLE.PhiBBPlusJet.combine_helper.combine_project import CombineProject
from DAZSLE.PhiBBPlusJet.combine_helper.region import Region
from ROOT import TFile, TGraph, TColor, gROOT, kOrange, kGreen, TCanvas, TLegend, TH1D, kBlack, kBlue
gROOT.SetBatch(True)
import math

# Load histograms
backgrounds = ["qcd", "tqq", "zqq", "wqq"]
signals = ["hqq125", "vbfhff125", "whqq125", "zhqq125"]

default_pt_categories = [[450., 500.], 	[500., 550.], [550., 600.], [600., 675.], [675., 800.], [800., 1000.]]

# Slice up the pt vs msd histogram
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
	return histogram_slices

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
		datacard_topdir = "/uscms/home/dryu/HBB/data/datacards/simple"
		datacard_dir = "{}/{}/{}".format(datacard_topdir, jet_type, signal)
		os.system("mkdir -pv {}".format(datacard_dir))
		project = CombineProject(datacard_dir)
		xvar = project.create_xvar("msd", [40., 600.])

		# Construct regions
		region_containers = {}
		for tag in ["vbf", "ggf"]:
			# Load 2D histograms
			histogram_file = TFile(histogram_files[tag], "READ")
			h2d = {}
			h2d_slices = {}
			for process in backgrounds + signals:
				h2d[process] = histogram_file.Get("{}_{}_{}".format(process, box, "pt_vs_msd"))
				if not h2d[process]:
					print "[simple_limits] ERROR : Couldn't load histogram {} from file {}".format("{}_{}_{}".format(process, box, "pt_vs_msd"), histogram_file.GetPath())
					sys.exit(1)

			# Normalize QCD to match data normalization
			#h2d_data_tmp = histogram_file.Get("data_obs_{}_pt_vs_msd".format(box))
			#int_data = h2d_data_tmp.Integral()
			#int_bkgds = {}
			#for process in backgrounds:
			#	int_bkgds[process] = h2d[process].Integral()
			#h2d["qcd"].Scale((int_data - sum([h2d[process].Integral() for process in backgrounds if process != "qcd"])) / h2d["qcd"].Integral())

			# Slice 2D histograms
			for process in backgrounds + signals:
				h2d_slices[process] = SliceSignalHistogram(h2d[process], pt_categories=default_pt_categories)

			# Data histogram = pseudodata from background
			for process in backgrounds:
				if not "data_obs" in h2d:
					h2d["data_obs"] = h2d[process].Clone()
					h2d["data_obs"].SetName("data_obs_{}_pt_vs_msd".format(box))
				else:
					h2d["data_obs"].Add(h2d[process])
			h2d_slices["data_obs"] = SliceSignalHistogram(h2d["data_obs"], pt_categories=default_pt_categories)

			for ptcat in xrange(len(default_pt_categories)):
				# Apply rho boundaries, and set data to zero if the prediction is zero
				if jet_type == "AK8":
					rho_range = [-6.0, -2.1]
				elif jet_type == "CA15":
					rho_range = [-4.7, -1.0]

				for xbin in xrange(1, h2d_slices["data_obs"][ptcat].GetNbinsX()+1):
					bin_mass = h2d_slices[process][ptcat].GetXaxis().GetBinCenter(xbin)
					bin_pt = default_pt_categories[ptcat][0] + (default_pt_categories[ptcat][1]-default_pt_categories[ptcat][0]) * 0.3
					bin_rho = 2. * math.log(bin_mass / bin_pt)
					if bin_rho < rho_range[0] or bin_rho > rho_range[1]:
						for process in backgrounds + ["data_obs", signal]:
							h2d_slices[process][ptcat].SetBinContent(xbin, 0.)
							h2d_slices[process][ptcat].SetBinError(xbin, 0.)

					# Background prediction in this bin
					bin_total_bkgd = 0.
					for process in backgrounds:
						bin_total_bkgd += h2d_slices[process][ptcat].GetBinContent(xbin)
					if bin_total_bkgd == 0:
						h2d_slices["data_obs"][ptcat].SetBinContent(xbin, 0)
						h2d_slices["data_obs"][ptcat].SetBinError(xbin, 0)

				region_name = "{}_ptcat{}".format(box, ptcat)
				region_containers[region_name] = Region(region_name)
				region_containers[region_name].add_xvar(xvar)

				for background in backgrounds:
					# Get pt vs msd histogram
					region_containers[region_name].add_background(background, h2d_slices[background][ptcat])
				region_containers[region_name].add_data("data_obs", h2d_slices["data_obs"][ptcat])
				region_containers[region_name].add_signal(signal, h2d_slices[signal][ptcat])
				project.add_region(region_containers[region_name])
		project.write()

		# Run all script
		with open("{}/run_all_limits.sh".format(datacard_topdir), 'w') as run_all_script:
			run_all_script.write("#!/bin/bash\n")
			for jet_type in ["AK8", "CA15"]:
				for signal in signals:
					datacard_dir = "{}/{}/{}".format(datacard_topdir, jet_type, signal)
					run_all_script.write("cd {}\n".format(datacard_dir))
					run_all_script.write("combine -M AsymptoticLimits datacard_total.txt >& limits.log\n")
					run_all_script.write("cd -\n".format(datacard_dir))
			run_all_script.close()

		# Sum of pass1 and pass2, i.e. no b-tag categories
		datacard_topdir = "/uscms/home/dryu/PhiBB2017/data/datacards/simple_nobtag"
		for jet_type in ["AK8", "CA15"]:
			for signal in signals:
				histogram_file = TFile("{}/histograms_SR_{}.root".format(histogram_dir, jet_type), "READ")
				datacard_dir = "{}/{}/{}".format(datacard_topdir, jet_type, signal)
				os.system("mkdir -pv {}".format(datacard_dir))
				project = CombineProject(datacard_dir)
				xvar = project.create_xvar("msd", [40., 600.])

				# Construct regions
				region_containers = {}

				# Load 2D histograms
				h2d = {}
				h2d_slices = {}
				for process in backgrounds + ["data_obs", signal]:
					h2d[process] = histogram_file.Get("{}_{}_{}".format(process, "pass1", "pt_vs_msd"))
					if not h2d[process]:
						print "[simple_limits] ERROR : Couldn't load histogram {} from file {}".format("{}_{}_{}".format(process, "pass1", "pt_vs_msd"), histogram_file.GetPath())
						sys.exit(1)
					h2d[process].Add(histogram_file.Get("{}_{}_{}".format(process, "pass2", "pt_vs_msd")))
					h2d_slices[process] = SliceSignalHistogram(h2d[process], pt_categories=default_pt_categories)

				for ptcat in xrange(len(default_pt_categories)):
					# Apply rho boundaries, and set data to zero if the prediction is zero
					if jet_type == "AK8":
						rho_range = [-6.0, -2.1]
					elif jet_type == "CA15":
						rho_range = [-4.7, -1.0]

					for xbin in xrange(1, h2d_slices["data_obs"][ptcat].GetNbinsX()+1):
						bin_mass = h2d_slices[process][ptcat].GetXaxis().GetBinCenter(xbin)
						bin_pt = default_pt_categories[ptcat][0] + (default_pt_categories[ptcat][1]-default_pt_categories[ptcat][0]) * 0.3
						bin_rho = 2. * math.log(bin_mass / bin_pt)
						if bin_rho < rho_range[0] or bin_rho > rho_range[1]:
							for process in backgrounds + ["data_obs", signal]:
								h2d_slices[process][ptcat].SetBinContent(xbin, 0.)
								h2d_slices[process][ptcat].SetBinError(xbin, 0.)

						# Background prediction in this bin
						bin_total_bkgd = 0.
						for process in backgrounds:
							bin_total_bkgd += h2d_slices[process][ptcat].GetBinContent(xbin)
						if bin_total_bkgd == 0:
							h2d_slices["data_obs"][ptcat].SetBinContent(xbin, 0)
							h2d_slices["data_obs"][ptcat].SetBinError(xbin, 0)

					region_name = "ptcat{}".format(ptcat)
					region_containers[region_name] = Region(region_name)
					region_containers[region_name].add_xvar(xvar)

					for background in backgrounds:
						# Get pt vs msd histogram
						region_containers[region_name].add_background(background, h2d_slices[background][ptcat])
					region_containers[region_name].add_data("data_obs", h2d_slices["data_obs"][ptcat])
					region_containers[region_name].add_signal(signal, h2d_slices[signal][ptcat])
					project.add_region(region_containers[region_name])
				project.write()

		# Run all script
		with open("{}/run_all_limits.sh".format(datacard_topdir), 'w') as run_all_script:
			run_all_script.write("#!/bin/bash\n")
			for jet_type in ["AK8", "CA15"]:
				for signal in signals:
					datacard_dir = "{}/{}/{}".format(datacard_topdir, jet_type, signal)
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

		for jet_type in ["AK8", "CA15"]:
			for signal_model in ["Sbb", "ZPrime"]:
				limit_whats = ["-2exp", "-1exp", "exp", "+1exp", "+2exp", "obs"]
				limits_btag = {}
				limits_nobtag = {}
				for limit_what in limit_whats:
					limits_btag[limit_what] = array.array('d', [])
					limits_nobtag[limit_what] = array.array('d', [])
				print "{} - {}".format(jet_type, signal_model)
				for mass in signal_model_masses[signal_model]:
					limit_file_btag = TFile("/uscms/home/dryu/PhiBB2017/data/datacards/simple/{}/{}{}/higgsCombineTest.AsymptoticLimits.mH120.root".format(jet_type, signal_model, int(mass)), "READ")
					print "Opening {}".format(limit_file_btag.GetPath())
					t_btag = limit_file_btag.Get("limit")
					for i, limit_what in enumerate(limit_whats):
						t_btag.GetEntry(i)
						limits_btag[limit_what].append(t_btag.GetLeaf("limit").GetValue(0))

					limit_file_nobtag = TFile("/uscms/home/dryu/PhiBB2017/data/datacards/simple_nobtag/{}/{}{}/higgsCombineTest.AsymptoticLimits.mH120.root".format(jet_type, signal_model, int(mass)), "READ")
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

				c = TCanvas("c_limitcomparison_{}_{}".format(signal_model, jet_type), "Limit comparison", 800, 600)
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

				c.SaveAs("/uscms/home/dryu/PhiBB2017/data/limits/figures/{}.pdf".format(c.GetName()))
			# End loop over signal models
		# End loop over jet types
