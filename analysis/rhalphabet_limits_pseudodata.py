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

pt_categories = [[450., 500.], 	[500., 550.], [550., 600.], [600., 675.], [675., 800.], [800., 1000.]]

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
			"SR":"/uscms/home/kkwok/work/Hbb/CMSSW_8_1_0/src/ZPrimePlusJet/fitting/PbbJet/ddb_Jan17/MC/data/hist_1DZbb_pt_scalesmear.root",
			"muonCR":"/uscms/home/kkwok/work/Hbb/CMSSW_8_1_0/src/ZPrimePlusJet/fitting/PbbJet/ddb_Jan17/MC/muonCR/hist_1DZbb_muonCR.root",
		}
		datacard_dir = "/uscms/home/dryu/HbbPlusJet/data/datacards/pseudodata"
		os.system("mkdir -pv {}".format(datacard_dir))
		project = CombineProject(datacard_dir)
		xvar = project.create_xvar("msd", [40., 201.])

		# Get histograms
		# Goal for maps: <category name, <process, histogram>>
		# category name = sr_ptcat#, or muoncr
		pass_histograms = {}
		fail_histograms = {}
		signal_regions = []

		# SR histograms
		histogram_file_SR = TFile(histogram_files["SR"], "READ")
		for box in ["pass", "fail"]:
			h2d = {}
			h2d_syst = {}

			# Load 2D histograms (ex.: ZPrime150_failn2_faildbtag_pt_vs_msd)
			for process in backgrounds + signals:
				hname = "{}_{}".format(process, box)
				if process in signals + ["zqq", "wqq"]:
					hname += "_matched"
				h2d[process] = histogram_file_SR.Get(hname)
				if not h2d[process]:
					print "[rhalphabet_limits] ERROR : Couldn't load histogram {} from file {}".format("{}_{}_{}".format(process, tag, "pt_vs_msd"), histogram_file_SR.GetPath())
					sys.exit(1)
				h2d[process].SetName("{}_SR_{}".format(process, box))
				h2d[process].SetDirectory(0)

				h2d_syst[process] = {}
				for syst in ["JER", "JES", "Pu", "trigger"]:
					for direction in ["Up", "Down"]:
						hname = "{}_{}_{}{}".format(process, box, syst, direction)
						h2d_syst[process][syst + direction] = histogram_file_SR.Get(hname)
						if not h2d_syst[process][syst + direction]:
							print "[rhalphabet_limits] ERROR : Couldn't load histogram {} from file {}".format(hname, histogram_file_SR.GetPath())
							sys.exit(1)
						h2d_syst[process][syst + direction].SetName("{}_SR_{}_{}{}".format(process, box, syst, direction))
						h2d_syst[process][syst + direction].SetDirectory(0)
			
			# Make pseudodata histogram
			for process in backgrounds + signals:
				if not "data_obs" in h2d:
					h2d["data_obs"] = h2d[process].Clone()
					h2d["data_obs"].SetName("data_obs_SR_{}".format(box))
					h2d["data_obs"].SetDirectory(0)
				else:
					h2d["data_obs"].Add(h2d[process])

			# Slice SR histograms
			histogram_slices = {}
			histogram_slices_syst = {}
			for process in backgrounds + signals + ["data_obs"]:
				histogram_slices[process] = SliceHistogram(h2d[process], pt_categories, rho_range=[-6.0, -2.1])
				histogram_slices_syst[process] = {}
				for syst in ["JER", "JES", "Pu", "trigger"]:
					for direction in ["Up", "Down"]:
						histogram_slices_syst[process][syst + direction] = SliceHistogram(h2d_syst[process][syst + direction], pt_categories, rho_range=[-6.0, -2.1])

			# Put the histogram slices into the dictionaries
			if box == "pass":
				hdict = pass_histograms
				hdict_syst = pass_histograms_syst
			else:
				hdict = fail_histograms
				hdict_syst = fail_histograms_syst
			for pt_index in len(pt_categories):
				region = "sr_ptcat{}".format(pt_index + 1)
				if not region in signal_regions:
					signal_regions.append(region)
				hdict[region] = {}
				hdict_syst[region] = {}
				for process in backgrounds + signals + ["data_obs"]:
					hdict[region][process] = histogram_slices[process][pt_index]
					hdict_syst[region][process] = {}
					for syst in ["JER", "JES", "Pu", "trigger"]:
						for direction in ["Up", "Down"]:
							hdict_syst["sr_ptcat{}".format(pt_index + 1)][process][syst + direction] = histogram_slices_syst[process][pt_index][syst + direction]
		histogram_file_SR.Close()

		# Get muon CR histograms
		histogram_file_muCR = TFile(histogram_files["muonCR"], "READ")
		pass_histograms["muonCR"] = {}
		fail_histograms["muonCR"] = {}
		pass_histograms_syst["muonCR"] = {}
		fail_histograms_syst["muonCR"] = {}
		for process in backgrounds + signals:
			hname_pass = "{}_pass".format(process) # Don't use matched histograms for muCR, because there's no rhalphabet QCD background estimation here
			pass_histograms["muonCR"][process] = histogram_file_muCR.Get(hname_pass)
			pass_histograms["muonCR"][process].SetDirectory(0)
			pass_histograms["muonCR"][process].SetName("{}_muonCR_pass".format(process))

			hname_fail = "{}_fail".format(process) # Don't use matched histograms for muCR, because there's no rhalphabet QCD background estimation here
			fail_histograms["muonCR"][process] = histogram_file_muCR.Get(hname_fail)
			fail_histograms["muonCR"][process].SetDirectory(0)
			fail_histograms["muonCR"][process].SetName("{}_muonCR_fail".format(process))

			pass_histograms_syst["muonCR"][process] = {}
			fail_histograms_syst["muonCR"][process] = {}
			for syst in ["JER", "JES", "Pu", "mutrigger", "muid"]:
				for direction in ["Up", "Down"]:
					hname_pass = "{}_pass_{}{}".format(process, syst, direction) # Don't use matched histograms for muCR, because there's no rhalphabet QCD background estimation here
					pass_histograms_syst["muonCR"][process] = histogram_file_muCR.Get(hname_pass)
					pass_histograms_syst["muonCR"][process].SetDirectory(0)
					pass_histograms_syst["muonCR"][process].SetName("{}_muonCR_pass".format(process))

					hname_fail = "{}_fail_{}{}".format(process, syst, direction) # Don't use matched histograms for muCR, because there's no rhalphabet QCD background estimation here
					fail_histograms_syst["muonCR"][process] = histogram_file_muCR.Get(hname_fail)
					fail_histograms_syst["muonCR"][process].SetDirectory(0)
					fail_histograms_syst["muonCR"][process].SetName("{}_muonCR_fail".format(process))

		# Again, make pseudodata
		for process in backgrounds + signals:
			if not "data_obs" in pass_histograms["muonCR"]:
				pass_histograms["muonCR"]["data_obs"] = pass_histograms["muonCR"][process].Clone()
				pass_histograms["muonCR"]["data_obs"].SetName("data_obs_muonCR_pass")
				fail_histograms["muonCR"]["data_obs"] = fail_histograms["muonCR"][process].Clone()
				fail_histograms["muonCR"]["data_obs"].SetName("data_obs_muonCR_fail")
			else:
				pass_histograms["muonCR"]["data_obs"].Add(pass_histograms["muonCR"][process])
				fail_histograms["muonCR"]["data_obs"].Add(fail_histograms["muonCR"][process])

		# For data histograms, if fail bin has 0 entries, set pass = 0
		for region, process_hist_dict in fail_histograms.iteritems():
			for process, hist in process_hist_dict.iteritems():
				for xbin in xrange(1, hist.GetNbinsX() + 1):
					if hist.GetBinContent(xbin) == 0:
						pass_histograms[region][process].SetBinContent(xbin, 0)

		
		#CreateCombineProject(region_names, pass_histograms, fail_histograms, jet_type=jet_type, datacard_directory=datacard_dir, signal_name=signal, data_name="data_obs", background_names=backgrounds, region_pts=region_pts)
		project = CombineProject(datacard_dir)	
		# Construct signal regions
		region_containers = {}
		for region_name in signal_regions:
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
			region_containers[region_name].add_rhbackground("qcd", y_value, n_x=args.nrho, n_y=args.npt, x_range=[-7., 0.], y_range=[400., 1000.])

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
		# Muon control region
		region_containers["muonCR"] = MuonControlRegion
		project.write()

		# Run all script
		with open("{}/run_limits.sh".format(datacard_dir), 'w') as run_all_script:
			run_all_script.write("#!/bin/bash\n")
			run_all_script.write("cd {}\n".format(datacard_dir))
			run_all_script.write("combine -M AsymptoticLimits datacard_total.txt >& limits.log\n")
			run_all_script.write("cd -\n".format(datacard_dir))
			run_all_script.close()


