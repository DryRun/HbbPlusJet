import os
import sys
from DAZSLE.DAZSLECommon.combine_helper.combine_project import CombineProject
from DAZSLE.DAZSLECommon.combine_helper.region import Region
from DAZSLE.DAZSLECommon.combine_helper.rhalphabet_region import RhalphabetRegion
from DAZSLE.DAZSLECommon.jet_mass_systematics import MassScaling, MassSmearing
from DAZSLE.DAZSLECommon.slice_histogram import SliceHistogram
from ROOT import TFile, TGraph, TColor, gROOT, kOrange, kGreen, TCanvas, TLegend, TH1D, kBlack, kBlue
gROOT.SetBatch(True)
import math

#import systematics

# Load histograms
backgrounds = ["qcd", "tqq", "zqq", "wqq"]
signals = ["hqq125", "tthqq125", "vbfhqq125", "whqq125", "zhqq125"]

default_pt_categories = [[450., 500.], 	[500., 550.], [550., 600.], [600., 675.], [675., 800.], [800., 1000.]]



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
		h2d_syst = {}
		h2d_slices_syst = {}
		for tag in ["ggf", "vbf"]:
			histogram_file = TFile(histogram_files[tag], "READ")
			for box in ["pass", "fail"]:
				region = "{}_{}".format(tag, box)
				h2d[region] = {}
				h2d_slices[region] = {}
				h2d_syst[region] = {}
				h2d_slices_syst[region] = {}

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

					h2d_syst[region][process] = {}
					h2d_slices_syst[region][process] = {}
					for syst in ["JER", "JES", "Pu", "trigger"]:
						for direction in ["Up", "Down"]:
							hname = "{}_{}_{}{}".format(process, box, syst, direction)
							h2d_syst[region][process][syst + direction] = histogram_file.Get(hname)
							if not h2d_syst[region][process][syst + direction]:
								print "[rhalphabet_limits] ERROR : Couldn't load histogram {} from file {}".format(hnamee, histogram_file.GetPath())
								sys.exit(1)
							h2d_syst[region][process][syst + direction].SetName("{}_{}_{}{}".format(process, region, syst, direction))
							h2d_syst[region][process][syst + direction].SetDirectory(0)

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
					h2d_slices[region][process] = SliceHistogram(h2d[region][process], pt_categories=default_pt_categories)
					for syst in ["JER", "JES", "Pu", "trigger"]:
						for direction in ["Up", "Down"]:
							h2d_slices_syst[region][process][syst + direction] = SliceHistogram(h2d_syst[region][process][syst + direction], pt_categories=default_pt_categories)

				# Data histogram = pseudodata from SM
				for process in backgrounds + signals:
					if not "data_obs" in h2d[region]:
						h2d[region]["data_obs"] = h2d[region][process].Clone()
						h2d[region]["data_obs"].SetName("data_obs_{}".format(region))
						h2d[region]["data_obs"].SetDirectory(0)
					else:
						h2d[region]["data_obs"].Add(h2d[region][process])
				h2d_slices[region]["data_obs"] = SliceHistogram(h2d[region]["data_obs"], pt_categories=default_pt_categories)
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
		pass_histograms_syst = {}
		fail_histograms_syst = {}
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
				pass_histograms_syst[region_name] = {}
				fail_histograms_syst[region_name] = {}
				for process in backgrounds + signals + ["data_obs"]:
					pass_histograms[region_name][process] = h2d_slices["{}_pass".format(tag)][process][ptcat-1]
					fail_histograms[region_name][process] = h2d_slices["{}_fail".format(tag)][process][ptcat-1]
					pass_histograms_syst[region_name][process] = {}
					fail_histograms_syst[region_name][process] = {}

				for process in backgrounds + signals:
					for syst in ["JER", "JES", "Pu", "trigger"]:
						for direction in ["Up", "Down"]:
							pass_histograms_syst[region_name][process][syst + direction] = h2d_slices_syst["{}_pass".format(tag)][process][syst + direction][ptcat-1]
							fail_histograms_syst[region_name][process][syst + direction] = h2d_slices_syst["{}_pass".format(tag)][process][syst + direction][ptcat-1]

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
			# def add_norm_systematic(self, syst_name, process, unc_value=None, unc_pass_fail=None):
			
			# Luminosity
			for process in backgrounds + signals:
				region_containers[region_name].add_norm_systematic("lumi", process, unc_value=1.025)

			# Higgs pT
			region_containers[region_name].add_norm_systematic("hqq125pt", "hqq125", unc_value=1.3)

			# N2 and double b-tag SFs
			for process in backgrounds + signals:
				if process in ["qcd", "tqq"]:
					continue
				region_containers[region_name].add_norm_systematic("veff", process, unc_value=1.043)

				# beff: need to derive fail region SF from pass+fail normalization constraint
				bbsf = 0.91
				dbbsf = 0.03
				if fail_histograms[region_name][process].Integral() > 0:
					dbbsf_fail = 1.0 - (dbbsf / bbsf) * pass_histograms[region_name][process].Integral() / fail_histograms[region_name][process].Integral()
				else:
					dbbsf_fail = None
				region_containers[region_name].add_norm_systematic("beff", process, unc_pass_fail=(1 + dbbsf / bbsf, dbbsf_fail))

			# Electroweak normalizations
			for process in ["wqq", "zqq"]:
				region_containers[region_name].add_norm_systematic("znormQ", process, unc_value=1.1)
				ptcat = region_cats[region_name]
				if ptcat in [1,2,3]:
					region_containers[region_name].add_norm_systematic("wznormEW", process, unc_value=1.05)
				elif ptcat in [4,5,6]:
					region_containers[region_name].add_norm_systematic("wznormEW", process, unc_value=1.15)
				znormEW_dict = {1:1.15, 2:1.15, 3:1.25, 4:1.35, 5:1.35, 6:1.35}
				region_containers[region_name].add_norm_systematic("znormEW", process, unc_value=znormEW_dict[ptcat])

			# JES, JER, Pu, trigger
			for systematic in ["JER", "JES", "Pu", "trigger"]:
				for process in backgrounds + signals:
					if process == "qcd":
						continue

					pass_norm = pass_histograms[region_name][process].Integral()
					pass_norm_up = pass_histograms_syst[region_name][process][syst + "Up"].Integral()
					pass_norm_down = pass_histograms_syst[region_name][process][syst + "Down"].Integral()
					if pass_norm > 0:
						pass_syst = 1.0 + (abs(pass_norm_up - pass_norm) + abs(pass_norm - pass_norm_down)) / (2. * pass_norm)
					else:
						pass_syst = "-"

					fail_norm = fail_histograms[region_name][process].Integral()
					fail_norm_up = fail_histograms_syst[region_name][process][syst + "Up"].Integral()
					fail_norm_down = fail_histograms_syst[region_name][process][syst + "Down"].Integral()
					if fail_norm > 0:
						fail_syst = 1.0 + (abs(fail_norm_up - fail_norm) + abs(fail_norm - fail_norm_down)) / (2. * fail_norm)
					else:
						fail_syst = "-"

					region_containers[region_name].add_norm_systematic(systematic, process, unc_pass_fail=(pass_syst, fail_syst))

			# Lepton veto
			for process in backgrounds + signals:
				if process == "qcd":
					continue
				region_containers[region_name].add_norm_systematic("muveto", process, unc_value=1.005)
				region_containers[region_name].add_norm_systematic("eleveto", process, unc_value=1.005)

			# Shape systematics 
			# Mass scale and resolution
			resonance_masses = {
				"wqq":80.385,
				"zqq":91.1876,
			}
			for signal in signals:
				resonance_masses[signal] = 125.
			for process in ["wqq", "zqq"] + signals:
				hists_massscale_pass = MassScaling(pass_histograms[region_name][process], resonance_masses[process])
				hists_massscale_fail = MassScaling(fail_histograms[region_name][process], resonance_masses[process])
				region_containers[region_name].add_shape_systematic("mass_scale", process, hists_massscale_pass, hists_massscale_fail, unc_value=0.1)

				hists_masssmear_pass = MassSmearing(pass_histograms[region_name][process])
				hists_masssmear_fail = MassSmearing(fail_histograms[region_name][process])
				region_containers[region_name].add_shape_systematic("mass_smear", process, hists_masssmear_pass, hists_masssmear_fail, unc_value=0.5)

			project.add_region(region_containers[region_name])
		project.write()

		# Run all script
		with open("{}/run_limits.sh".format(datacard_dir), 'w') as run_all_script:
			run_all_script.write("#!/bin/bash\n")
			run_all_script.write("cd {}\n".format(datacard_dir))
			run_all_script.write("combine -M AsymptoticLimits datacard_total.txt >& limits.log\n")
			run_all_script.write("cd -\n".format(datacard_dir))
			run_all_script.close()


