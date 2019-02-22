import os
import sys
import ROOT
from DAZSLE.PhiBBPlusJet.analysis_base import AnalysisBase
import DAZSLE.PhiBBPlusJet.analysis_configuration as config
import DAZSLE.PhiBBPlusJet.event_selections as event_selections
from DAZSLE.PhiBBPlusJet.event_selector import *
from DAZSLE.ZPrimePlusJet.xbb_config import analysis_parameters as params
from math import ceil, sqrt,floor
import array

import ROOT
from ROOT import *
gInterpreter.Declare("#include \"MyTools/RootUtils/interface/SeabornInterface.h\"")
gInterpreter.Declare("#include \"MyTools/RootUtils/interface/HistogramManager.h\"")
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libMyToolsRootUtils.so"))
gInterpreter.Declare("#include \"MyTools/AnalysisTools/interface/EventSelector.h\"")
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libMyToolsAnalysisTools.so"))
ROOT.gInterpreter.Declare("#include \"DAZSLE/PhiBBPlusJet/interface/BaconData.h\"")
ROOT.gInterpreter.Declare("#include \"DAZSLE/PhiBBPlusJet/interface/BaconEventCutFunctions.h\"")
gSystem.Load(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libDAZSLEPhiBBPlusJet.so"))
#from ROOT import gInterpreter, gSystem, gROOT, gStyle, Root, TCanvas, TLegend, TH1F, TFile, TGraphErrors
gROOT.SetBatch(ROOT.kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
seaborn = Root.SeabornInterface()
seaborn.Initialize()

# Enums
#from BaconData import kAK8, kCA15, kPt, kDbtag, kN2DDT

class Histograms(AnalysisBase):
	def __init__(self, sample_name, tree_name="otree", jet_type="AK8", jet_ordering="pt"):
		super(Histograms, self).__init__(tree_name=tree_name)
		self._data = BaconData(self._chain)
		self._output_path = ""
		self._sample_name = sample_name
		self._input_nevents = 0
		self._processed_events = 0
		self._pt_bins = array.array("d", [450., 500.,550.,600.,675.,800.,1000.])
		if jet_ordering == "pt":
			self._jet_ordering = BaconData.kPt
		elif jet_ordering == "dbtag":
			self._jet_ordering = BaconData.kDbtag
		elif jet_ordering == "n2ddt":
			self._jet_ordering = BaconData.kN2DDT
		self._jet_type = jet_type
		if self._jet_type == "AK8":
			jet_type_enum = BaconData.kAK8
		elif self._jet_type == "CA15":
			jet_type_enum = BaconData.kCA15
		self._data.SetJetSelection(jet_type_enum, self._jet_ordering)

		# Cuts
		self._n2ddt_wps = [0.05, 0.15, 0.26]
		self._dbtag_cuts= [0.7, 0.8, 0.9]
		#self._dbtagcut = 0.8
		#self._dbtagcut_loose = 0.7
		self._dcsv_min = -999.
		self._data_source = "data"
		self._prescale = -1

		# Systematics
		# Weight systematics: these only affect the weights used to fill histograms, so can easily be filled in normal running
		self._weight_systematics = {
			"SR":["TriggerUp", "TriggerDown", "PUUp", "PUDown"],
			"Preselection":["TriggerUp", "TriggerDown", "PUUp", "PUDown"],
			"muCR":["MuTriggerUp", "MuTriggerDown", "MuIDUp", "MuIDDown", "MuIsoUp", "MuIsoDown", "PUUp", "PUDown"]
		}
		self._weight_systematics["SR_matched"] = self._weight_systematics["SR"]
		# Jet systematics: these affect the jet pT, so modify the event selection
		self._jet_systematics = ["JESUp", "JESDown", "JERUp", "JERDown"]

	def set_prescale(self, prescale):
		self._prescale = prescale

	def set_data_source(self, data_source):
		print "[Histograms::set_data_source] INFO : Setting data source to " + data_source
		if not data_source in ["data", "simulation"]:
			print "[Histograms] ERROR : Data source must be data or simulation."
			sys.exit(1)
		self._data_source = data_source

	# Overload add_file to extract the number of input events to the skims, stored in histogram NEvents in the same file as the trees
	def add_file(self, filename):
		super(Histograms, self).add_file(filename)
		f = ROOT.TFile.Open(filename, "READ")
		if f.Get("NEvents").Integral() == 0:
			print "[Histograms::add_file] ERROR : NEvents.Integral() == 0 for file " + filename
			sys.exit(1)
		self._input_nevents += f.Get("NEvents").Integral()
		f.Close()

	def set_output_path(self, output_path):
		self._output_path = output_path
		os.system("mkdir -pv {}".format(os.path.dirname(self._output_path)))

	def initialize_weight_tools(self):
		# Pileup weight stuff
		f_pu = TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/puWeights_All.root", "read")
		self._h_pu_weight = f_pu.Get("puw")
		self._h_pu_weight.SetDirectory(0)
		self._h_pu_weight_up = f_pu.Get("puw_p")
		self._h_pu_weight_up.SetDirectory(0)
		self._h_pu_weight_down = f_pu.Get("puw_m")
		self._h_pu_weight_down.SetDirectory(0)
		f_pu.Close()

		# Trigger efficiency weight stuff
		if self._jet_type == "AK8":
			f_trig = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/RUNTriggerEfficiencies_AK8_SingleMuon_Run2016_V2p1_v03.root", "read")
			self._trig_den = f_trig.Get("DijetTriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtDenom_cutJet")
			self._trig_num = f_trig.Get("DijetTriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtPassing_cutJet")
		elif self._jet_type == "CA15":
			f_trig = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/RUNTriggerEfficiencies_CA15_SingleMuon_Run2016_V2p4_v08.root", "read")
			self._trig_den = f_trig.Get("DijetCA15TriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtDenom_cutJet")
			self._trig_num = f_trig.Get("DijetCA15TriggerEfficiencySeveralTriggers/jet1SoftDropMassjet1PtPassing_cutJet")
		self._trig_den.SetDirectory(0)
		self._trig_num.SetDirectory(0)
		self._trig_den.RebinX(2)
		self._trig_num.RebinX(2)
		self._trig_den.RebinY(5)
		self._trig_num.RebinY(5)
		self._trig_eff = ROOT.TEfficiency()
		if (ROOT.TEfficiency.CheckConsistency(self._trig_num, self._trig_den)):
			self._trig_eff = ROOT.TEfficiency(self._trig_num, self._trig_den)
			self._trig_eff.SetDirectory(0)
		f_trig.Close()

		# get muon trigger efficiency object

		lumi_GH = 16.146
		lumi_BCDEF = 19.721
		lumi_total = lumi_GH + lumi_BCDEF

		f_mutrig_GH = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/EfficienciesAndSF_Period4.root", "read")
		self._mutrig_eff_GH = f_mutrig_GH.Get("Mu50_OR_TkMu50_PtEtaBins/efficienciesDATA/pt_abseta_DATA")
		self._mutrig_eff_GH.Sumw2()
		self._mutrig_eff_GH.SetDirectory(0)
		f_mutrig_GH.Close()

		f_mutrig_BCDEF = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/EfficienciesAndSF_RunBtoF.root", "read")
		self._mutrig_eff_BCDEF = f_mutrig_BCDEF.Get("Mu50_OR_TkMu50_PtEtaBins/efficienciesDATA/pt_abseta_DATA")
		self._mutrig_eff_BCDEF.Sumw2()
		self._mutrig_eff_BCDEF.SetDirectory(0)
		f_mutrig_BCDEF.Close()

		self._mutrig_eff = self._mutrig_eff_GH.Clone('pt_abseta_DATA_mutrig_ave')
		self._mutrig_eff.Scale(lumi_GH / lumi_total)
		self._mutrig_eff.Add(self._mutrig_eff_BCDEF, lumi_BCDEF / lumi_total)

		# get muon ID efficiency object

		f_muid_GH = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/EfficienciesAndSF_GH.root", "read")
		self._muid_eff_GH = f_muid_GH.Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/pt_abseta_DATA")
		self._muid_eff_GH.Sumw2()
		self._muid_eff_GH.SetDirectory(0)
		f_muid_GH.Close()

		f_muid_BCDEF = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/EfficienciesAndSF_BCDEF.root", "read")
		self._muid_eff_BCDEF = f_muid_BCDEF.Get(
			"MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/pt_abseta_DATA")
		self._muid_eff_BCDEF.Sumw2()
		self._muid_eff_BCDEF.SetDirectory(0)
		f_muid_BCDEF.Close()

		self._muid_eff = self._muid_eff_GH.Clone('pt_abseta_DATA_muid_ave')
		self._muid_eff.Scale(lumi_GH / lumi_total)
		self._muid_eff.Add(self._muid_eff_BCDEF, lumi_BCDEF / lumi_total)

		# get muon ISO efficiency object

		f_muiso_GH = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/EfficienciesAndSF_ISO_GH.root", "read")
		self._muiso_eff_GH = f_muiso_GH.Get("LooseISO_LooseID_pt_eta/efficienciesDATA/pt_abseta_DATA")
		self._muiso_eff_GH.Sumw2()
		self._muiso_eff_GH.SetDirectory(0)
		f_muiso_GH.Close()

		f_muiso_BCDEF = ROOT.TFile.Open("$CMSSW_BASE/src/DAZSLE/ZPrimePlusJet/analysis/ggH/EfficienciesAndSF_ISO_BCDEF.root", "read")
		self._muiso_eff_BCDEF = f_muiso_BCDEF.Get("LooseISO_LooseID_pt_eta/efficienciesDATA/pt_abseta_DATA")
		self._muiso_eff_BCDEF.Sumw2()
		self._muiso_eff_BCDEF.SetDirectory(0)
		f_muiso_BCDEF.Close()

		self._muiso_eff = self._muiso_eff_GH.Clone('pt_abseta_DATA_muiso_ave')
		self._muiso_eff.Scale(lumi_GH / lumi_total)
		self._muiso_eff.Add(self._muiso_eff_BCDEF, lumi_BCDEF / lumi_total)

	def start(self):
		# Event selections
		self._selections = []
		self._event_selectors = {}

		### Preselection ###
		self._selections.append("Preselection")
		self._event_selectors["Preselection"] = EventSelector("Preselection")

		@add_cut(self._event_selectors["Preselection"])
		@add_nm1_hist(self._event_selectors["Preselection"], "pt", "p_{T} [GeV]", 120, 0., 1200.)
		def min_pt(self, event):
			self._return_data["min_pt"]["pt"] = event.SelectedJet_pt
			return event.SelectedJet_pt > 400.

		@add_cut(self._event_selectors["Preselection"])
		@add_nm1_hist(self._event_selectors["Preselection"], "msd_puppi", "p_{T} [GeV]", 120, 0., 1200.)
		def min_msd(self, event):
			self._return_data["min_msd"]["msd_puppi"] = event.SelectedJet_msd_puppi
			return event.SelectedJet_msd_puppi > 40.

		### Signal region ###
		self._selections.append("SR")
		self._event_selectors["SR"] = EventSelector("SR")

		@add_cut(self._event_selectors["SR"])
		def jetID(self, event):
			return event.SelectedJet_isTightVJet == 1

		@add_cut(self._event_selectors["SR"])
		@add_nm1_hist(self._event_selectors["SR"], "pt", "p_{T} [GeV]", 120, 0., 1200.)
		def min_pt(self, event):
			self._return_data["min_pt"]["pt"] = event.SelectedJet_pt
			return event.SelectedJet_pt > 450.

		@add_cut(self._event_selectors["SR"])
		@add_nm1_hist(self._event_selectors["SR"], "msd", "m_{SD}^{PUPPI} [GeV]", 80, 40., 600.)
		def min_msd(self, event):
			self._return_data["min_msd"]["msd"] = event.SelectedJet_msd_puppi
			return event.SelectedJet_msd_puppi > 40.

		@add_cut(self._event_selectors["SR"])
		@add_nm1_hist(self._event_selectors["SR"], "neleLoose", "N_{e}", 11, -0.5, 10.5)
		def electronveto(self, event):
			self._return_data["electronveto"]["neleLoose"] = event.neleLoose
			return event.neleLoose == 0

		@add_cut(self._event_selectors["SR"])
		@add_nm1_hist(self._event_selectors["SR"], "nmuLoose", "N_{#mu}", 11, -0.5, 10.5)
		def muonveto(self, event):
			self._return_data["muonveto"]["nmuLoose"] = event.nmuLoose
			return event.nmuLoose == 0

		@add_cut(self._event_selectors["SR"])
		@add_nm1_hist(self._event_selectors["SR"], "ntau", "N_{#tau}", 11, -0.5, 10.5)
		def tauveto(self, event):
			self._return_data["tauveto"]["ntau"] = event.ntau
			return event.ntau==0

		@add_cut(self._event_selectors["SR"])
		@add_nm1_hist(self._event_selectors["SR"], "pfmet", "pfmet [GeV]", 25, 0., 500.)
		def max_pfmet(self, event):
			self._return_data["max_pfmet"]["pfmet"] = event.pfmet
			return event.pfmet<140.

		### Signal region, with additional truth-reco V matching ###
		if self._data_source == "simulation":
			self._selections.append("SR_matched")
			self._event_selectors["SR_matched"] = EventSelector("SR_matched")

			@add_cut(self._event_selectors["SR_matched"])
			@add_nm1_hist(self._event_selectors["SR_matched"], "Vmatch_dphi", "#Delta #phi(V_{truth}, V_{rec})", 140, -7., 7.)
			@add_nm1_hist(self._event_selectors["SR_matched"], "Vmatch_dR", "#Delta R(V_{truth}, V_{rec})", 100, 0., 10.)
			@add_nm1_hist(self._event_selectors["SR_matched"], "Vmatch_dpt", "(p_{T}^{gen}-p_{T}^{rec})/p_{T}^{gen}", 200, -10., 10.)
			def Vmatch(self, event):
				if event.genVPt > 0 and event.genVMass > 0:
					matching_dphi = abs(math.acos(math.cos(event.genVPhi - event.SelectedJet_phi)))
					matching_deta = abs(event.genVEta - event.SelectedJet_eta)
					matching_dpt = abs(event.genVPt - event.SelectedJet_pt) / event.genVPt
					matching_dmass = abs(event.genVMass - event.SelectedJet_msd_puppi) / event.genVMass
					vmatched = matching_dphi < 0.8 and matching_dpt < 0.5 and matching_dmass < 0.3
					self._return_data["Vmatch"]["Vmatch_dphi"] = matching_dphi
					self._return_data["Vmatch"]["Vmatch_dR"] = (matching_deta**2 + matching_dphi**2)**0.5
					self._return_data["Vmatch"]["Vmatch_dpt"] = (event.genVPt - event.SelectedJet_pt) / event.genVPt

					return (event.genVPt>0 and event.genVMass>0 and matching_dphi < 0.8 and matching_dpt < 0.5 and matching_dmass < 0.3)
				else:
					self._return_data["Vmatch"]["Vmatch_dphi"] = -1.e20
					self._return_data["Vmatch"]["Vmatch_dR"] = -1.e20
					self._return_data["Vmatch"]["Vmatch_dpt"] = -1.e20
					return False

			@add_cut(self._event_selectors["SR_matched"])
			def jetID(self, event):
				return event.SelectedJet_isTightVJet == 1

			@add_cut(self._event_selectors["SR_matched"])
			@add_nm1_hist(self._event_selectors["SR_matched"], "pt", "p_{T} [GeV]", 120, 0., 1200.)
			def min_pt(self, event):
				self._return_data["min_pt"]["pt"] = event.SelectedJet_pt
				return event.SelectedJet_pt > 450.

			@add_cut(self._event_selectors["SR_matched"])
			@add_nm1_hist(self._event_selectors["SR_matched"], "msd", "m_{SD}^{PUPPI} [GeV]", 80, 40., 600.)
			def min_msd(self, event):
				self._return_data["min_msd"]["msd"] = event.SelectedJet_msd_puppi
				return event.SelectedJet_msd_puppi > 40.

			@add_cut(self._event_selectors["SR_matched"])
			@add_nm1_hist(self._event_selectors["SR_matched"], "neleLoose", "N_{e}", 11, -0.5, 10.5)
			def electronveto(self, event):
				self._return_data["electronveto"]["neleLoose"] = event.neleLoose
				return event.neleLoose == 0

			@add_cut(self._event_selectors["SR_matched"])
			@add_nm1_hist(self._event_selectors["SR_matched"], "nmuLoose", "N_{#mu}", 11, -0.5, 10.5)
			def muonveto(self, event):
				self._return_data["muonveto"]["nmuLoose"] = event.nmuLoose
				return event.nmuLoose == 0

			@add_cut(self._event_selectors["SR_matched"])
			@add_nm1_hist(self._event_selectors["SR_matched"], "ntau", "N_{#tau}", 11, -0.5, 10.5)
			def tauveto(self, event):
				self._return_data["tauveto"]["ntau"] = event.ntau
				return event.ntau==0

			@add_cut(self._event_selectors["SR_matched"])
			@add_nm1_hist(self._event_selectors["SR_matched"], "pfmet", "pfmet [GeV]", 25, 0., 500.)
			def max_pfmet(self, event):
				self._return_data["max_pfmet"]["pfmet"] = event.pfmet
				return event.pfmet<140.

		for systematic in self._jet_systematics:
			selector_syst_name = "SR_{}".format(systematic)
			self._selections.append(selector_syst_name)
			self._event_selectors[selector_syst_name] =  EventSelector(selector_syst_name)
			@add_cut(self._event_selectors[selector_syst_name])
			def jetID(self, event):
				return (event.SelectedJet_isTightVJet == 1)
			@add_cut(self._event_selectors[selector_syst_name])
			def min_pt(self, event):
				return (event.SelectedJet_pt > 450.)
			@add_cut(self._event_selectors[selector_syst_name])
			def min_msd(self, event):
				return (event.SelectedJet_msd_puppi > 40.)
			@add_cut(self._event_selectors[selector_syst_name])
			def electronveto(self, event):
				return (event.neleLoose==0)
			@add_cut(self._event_selectors[selector_syst_name])
			def muonveto(self, event):
				return (event.nmuLoose==0)
			@add_cut(self._event_selectors[selector_syst_name])
			def tauveto(self, event):
				return (event.ntau==0)
			@add_cut(self._event_selectors[selector_syst_name])
			def max_pfmet(self, event):
				return (event.pfmet<140.)

			if self._data_source == "simulation":
				selector_syst_name = "SR_{}_matched".format(systematic)
				self._selections.append(selector_syst_name)
				self._event_selectors[selector_syst_name] = EventSelector(selector_syst_name)
				
				@add_cut(self._event_selectors[selector_syst_name])
				def Vmatch(self, event):
					if event.genVPt > 0 and event.genVMass > 0:
						matching_dphi = abs(math.acos(math.cos(event.genVPhi - event.SelectedJet_phi))) 
						matching_dpt = abs(event.genVPt - event.SelectedJet_pt) / event.genVPt 
						matching_dmass = abs(event.genVMass - event.SelectedJet_msd_puppi) / event.genVMass 
						vmatched = matching_dphi < 0.8 and matching_dpt < 0.5 and matching_dmass < 0.3 
						return (event.genVPt>0 and event.genVMass>0 and matching_dphi < 0.8 and matching_dpt < 0.5 and matching_dmass < 0.3)
					else:
						return False

				@add_cut(self._event_selectors[selector_syst_name])
				def jetID(self, event):
				 	return (event.SelectedJet_isTightVJet == 1)

				@add_cut(self._event_selectors[selector_syst_name])
				def min_pt(self, event):
				 	return (event.SelectedJet_pt > 450.)

				@add_cut(self._event_selectors[selector_syst_name])
				def min_msd(self, event):
				 	return (event.SelectedJet_msd_puppi > 40.)

				@add_cut(self._event_selectors[selector_syst_name])
				def electronveto(self, event):
					return (event.neleLoose==0)

				@add_cut(self._event_selectors[selector_syst_name])
				def muonveto(self, event):
				 	return (event.nmuLoose==0)

				@add_cut(self._event_selectors[selector_syst_name])
				def tauveto(self, event):
				 	return (event.ntau==0)

				@add_cut(self._event_selectors[selector_syst_name])
				def max_pfmet(self, event):
				 	return (event.pfmet<140.)


		# Histograms
		self._histograms = ROOT.Root.HistogramManager()
		self._histograms.AddPrefix("h_")
		self._histograms.AddTH1F("input_nevents", "input_nevents", "", 1, -0.5, 0.5)
		self._histograms.GetTH1F("input_nevents").SetBinContent(1, self._input_nevents)
		self._histograms.AddTH1D("processed_nevents", "processed_nevents", "", 1, -0.5, 0.5)
		self._histograms.AddTH1D("inclusive_pt", "inclusive_pt", "p_{T} [GeV]", 100, 0., 1000.)
		self._histograms.AddTH1D("inclusive_eta", "inclusive_eta", "#eta", 100, -5., 5.)
		self._histograms.AddTH1D("inclusive_msd", "inclusive_msd", "", 100, 0., 1000.)
		self._histograms.AddTH1D("inclusive_n2ddt", "inclusive_n2ddt", "", 40, -2., 2.)
		self._histograms.AddTH1D("inclusive_dcsv", "inclusive_dcsv", "", 200, -1., 1.)

		# Histograms for each event selection
		self._selection_histograms = {}
		self._boxes = ["all"]
		self._wps = [] # List of (n2ddt wp, dbtag cut)
		# "pass1", "pass2", "fail1", "fail2", "pass1loose", "pass2loose"]
		for n2ddt_wp in self._n2ddt_wps:
			for dbtag_cut in self._dbtag_cuts:
				self._wps.append((n2ddt_wp, dbtag_cut))
				box_suffix = "_n2wp{}_dbtag{}".format(n2ddt_wp, dbtag_cut)
				self._boxes.append("passn2_passdbtag{}".format(box_suffix))
				self._boxes.append("passn2_faildbtag{}".format(box_suffix))
				self._boxes.append("failn2_passdbtag{}".format(box_suffix))
				self._boxes.append("failn2_faildbtag{}".format(box_suffix))

		for selection in self._selections:
			self._selection_histograms[selection] = ROOT.Root.HistogramManager()
			self._selection_histograms[selection].AddPrefix("h_{}_".format(selection))

			for box in self._boxes:
				self._selection_histograms[selection].AddTH2D("{}_pt_vs_msd".format(box), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", 240, 0., 1200.)
				self._selection_histograms[selection].AddTH2D("{}_pt_vs_msd_unweighted".format(box), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", 240, 0., 1200.)
				self._selection_histograms[selection].AddTH1D("{}_nevents".format(box), "{}_nevents".format(box), "", 1, -0.5, 0.5)
				self._selection_histograms[selection].AddTH1D("{}_nevents_weighted".format(box), "{}_nevents_weighted".format(box), "", 1, -0.5, 0.5)
				self._selection_histograms[selection].AddTH1D("{}_pfmet".format(box), "PF MET", "PF MET [GeV]", 200, 0., 1000.)
				self._selection_histograms[selection].AddTH1D("{}_dcsv".format(box), "dcsv", "dcsv", 200, -1., 1.)
				self._selection_histograms[selection].AddTH1D("{}_n2ddt".format(box), "n2ddt", "N_{2}^{DDT}", 160, -1.0, 1.0)
				self._selection_histograms[selection].AddTH1D("{}_n2".format(box), "n2", "N_{2}", 160, -1.0, 1.0)
				self._selection_histograms[selection].AddTH1D("{}_tau21ddt".format(box), "tau21ddt", "tau21ddt", 20, -0.5, 0.5)
				self._selection_histograms[selection].AddTH1D("{}_pt".format(box), "pt", "pt", 400, 0., 2000.)
				self._selection_histograms[selection].AddTH1D("{}_msd".format(box), "msd", "msd", 85, 5, 600)
				self._selection_histograms[selection].AddTH1D("{}_eta".format(box), "eta", "eta", 60, -3., 3.)
				self._selection_histograms[selection].AddTH1D("{}_rho".format(box), "rho", "rho", 80, -8., 0.)
				self._selection_histograms[selection].AddTH2D("{}_met_vs_msd".format(box), "met_msd", "m_{SD} [GeV]", 40, 40, 600, "E_{T}^{miss} [GeV]", 25, 0., 500.)
				if selection in self._weight_systematics:
					for systematic in self._weight_systematics[selection]:
						self._selection_histograms[selection].AddTH2D("{}_{}_pt_vs_msd".format(systematic, box), "; {} m_{{SD}}^{{PUPPI}} (GeV); {} p_{{T}} (GeV)".format(self._jet_type, self._jet_type), "m_{SD}^{PUPPI} [GeV]", 80, 40, 600, "p_{T} [GeV]", 240, 0., 1200.)

			self._selection_histograms[selection].AddTH3D("dbtag_vs_pt_vs_msd", "dbtag_vs_pt_vs_msd", 
				"m_{SD} [GeV]", 80, 40, 600,
				"p_{T} [GeV]", 12, 400, 1000,
				"Double-b", 110, -1.1, 1.1)
			self._selection_histograms[selection].AddTH3D("n2ddt_vs_pt_vs_msd", "n2ddt_vs_pt_vs_msd", 
				"m_{SD} [GeV]", 40, 40, 600,
				"p_{T} [GeV]", 12, 400, 1000,
				"N_{2}^{DDT}", 40, -0.5, 0.5)
			self._selection_histograms[selection].AddTH2D("nparticles_vs_n2", "N_{particles} vs N_{2}^{1}", "N_{2}^{1}", 40, -1., 1., "N_{particles}", 41, -0.5, 40.5)
		self.initialize_weight_tools()


	def run(self, max_nevents=-1, first_event=0):
		if max_nevents > 0:
			limit_nevents = min(max_nevents, self._chain.GetEntries())
		else:
			limit_nevents = self._chain.GetEntries()

		n_checkpoints = 20
		print_every = int(ceil(1. * limit_nevents / n_checkpoints))

		print "[EventSelectionHistograms::run] INFO : Running loop over tree from event {} to {}".format(first_event, limit_nevents - 1)

		self.start_timer()
		for entry in xrange(first_event, limit_nevents):
			self.print_progress(entry, first_event, limit_nevents, print_every)
			self._data.GetEntry(entry)

			# Prescale before anything
			if self._prescale > 0:
				if self._data.evtNum % self._prescale != 0:
					continue

			self._histograms.GetTH1D("processed_nevents").Fill(0)
			self._processed_events += 1

			npu = min(self._data.npu, 49.5)
			pu_weight = self._h_pu_weight.GetBinContent(self._h_pu_weight.FindBin(npu))
			pu_weight_up = self._h_pu_weight_up.GetBinContent(self._h_pu_weight_up.FindBin(npu))
			pu_weight_down = self._h_pu_weight_down.GetBinContent(self._h_pu_weight_down.FindBin(npu))

			k_vjets = 1.
			k_ttbar = 1.
			w_scale = {
				(0, 500):1.0,
				(500, 600):1.0,
				(600, 700):1.0,
				(700, 800):1.2,
				(800, 900):1.25,
				(900, 1000):1.25,
				(1000, 3000):1.0
			}
			if self._sample_name == 'wqq' or self._sample_name == 'W':
				k_vjets = self._data.kfactor * 1.35  # ==1 for not V+jets events
				for pt_range, w_sf in w_scale.iteritems():
					if pt_range[0] < self._data.genVPt < pt_range[1]:
						k_vjets *= w_sf
			elif self._sample_name == 'zqq' or self._sample_name == 'DY':
				k_vjets = self._data.kfactor * 1.45  # ==1 for not V+jets events
			elif self._sample_name == 'tqq':
				k_ttbar = self._data.topPtWeight

			for selection in self._selections:
				# Get weights
				if self._data_source == "data":
					event_weight = 1.
					event_weight_syst = {}
					if "SR" in selection or "Preselection" in selection or "N2CR" in selection:
						event_weight_syst["TriggerUp"] = 1.
						event_weight_syst["TriggerDown"] = 1.
						event_weight_syst["PUUp"] = 1.
						event_weight_syst["PUDown"] = 1.
					elif "muCR" in selection:
						event_weight_syst["MuTriggerUp"] = 1.
						event_weight_syst["MuTriggerDown"] = 1.
						event_weight_syst["MuIDUp"] = 1.
						event_weight_syst["MuIDDown"] = 1.
						event_weight_syst["MuIsoUp"] = 1.
						event_weight_syst["MuIsoDown"] = 1.
						event_weight_syst["PUUp"] = 1.
						event_weight_syst["PUDown"] = 1.
					else:
						print "[event_selection_histograms::run] ERROR : selection {} is not known. Fix!".format(selection)
						sys.exit(1)
				else:
					if "SR" in selection or "Preselection" in selection or "N2CR" in selection:
						if self._jet_type == "AK8":
							trigger_mass = min(self._data.AK8Puppijet0_msd, 300.)
							trigger_pt = max(200., min(self._data.AK8Puppijet0_pt, 1000.))
						elif self._jet_type == "CA15":
							trigger_mass = min(self._data.CA15Puppijet0_msd, 300.)
							trigger_pt = max(200., min(self._data.CA15Puppijet0_pt, 1000.))
						trigger_weight = self._trig_eff.GetEfficiency(self._trig_eff.FindFixBin(trigger_mass, trigger_pt))
						trigger_weight_up = trigger_weight + self._trig_eff.GetEfficiencyErrorUp(self._trig_eff.FindFixBin(trigger_mass, trigger_pt))
						trigger_weight_down = trigger_weight - self._trig_eff.GetEfficiencyErrorLow(
							self._trig_eff.FindFixBin(trigger_mass, trigger_pt))
						if trigger_weight <= 0 or trigger_weight_down <= 0 or trigger_weight_up <= 0:
							#print 'trigger_weights are %f, %f, %f, setting all to 1' % (trigger_weight, trigger_weight_up, trigger_weight_down)
							trigger_weight = 1
							trigger_weight_down = 1
							trigger_weight_up = 1

						event_weight = pu_weight * k_vjets * trigger_weight * k_ttbar
						event_weight_syst = {}
						event_weight_syst["TriggerUp"] = pu_weight * k_vjets * k_ttbar *trigger_weight_up
						event_weight_syst["TriggerDown"] = pu_weight * k_vjets * k_ttbar * trigger_weight_down
						event_weight_syst["PUUp"] = pu_weight_up * k_vjets * k_ttbar * trigger_weight
						event_weight_syst["PUDown"] = pu_weight_down * k_vjets * k_ttbar * trigger_weight

					elif "muCR" in selection:
						mutrigweight = 1
						mutrigweightDown = 1
						mutrigweightUp = 1
						if self._data.nmuLoose > 0:
							muPtForTrig = max(52., min(self._data.vmuoLoose0_pt, 700.))
							muEtaForTrig = min(abs(self._data.vmuoLoose0_eta), 2.3)
							mutrigweight = self._mutrig_eff.GetBinContent(self._mutrig_eff.FindBin(muPtForTrig, muEtaForTrig))
							mutrigweightUp = mutrigweight + self._mutrig_eff.GetBinError(
								self._mutrig_eff.FindBin(muPtForTrig, muEtaForTrig))
							mutrigweightDown = mutrigweight - self._mutrig_eff.GetBinError(
								self._mutrig_eff.FindBin(muPtForTrig, muEtaForTrig))
							if mutrigweight <= 0 or mutrigweightDown <= 0 or mutrigweightUp <= 0:
								print 'mutrigweights are %f, %f, %f, setting all to 1' % (
								mutrigweight, mutrigweightUp, mutrigweightDown)
								mutrigweight = 1
								mutrigweightDown = 1
								mutrigweightUp = 1

						muidweight = 1
						muidweightDown = 1
						muidweightUp = 1
						if self._data.nmuLoose > 0:
							muPtForId = max(20., min(self._data.vmuoLoose0_pt, 100.))
							muEtaForId = min(abs(self._data.vmuoLoose0_eta), 2.3)
							muidweight = self._muid_eff.GetBinContent(self._muid_eff.FindBin(muPtForId, muEtaForId))
							muidweightUp = muidweight + self._muid_eff.GetBinError(self._muid_eff.FindBin(muPtForId, muEtaForId))
							muidweightDown = muidweight - self._muid_eff.GetBinError(self._muid_eff.FindBin(muPtForId, muEtaForId))
							if muidweight <= 0 or muidweightDown <= 0 or muidweightUp <= 0:
								print 'muidweights are %f, %f, %f, setting all to 1' % (muidweight, muidweightUp, muidweightDown)
								muidweight = 1
								muidweightDown = 1
								muidweightUp = 1

						muisoweight = 1
						muisoweightDown = 1
						muisoweightUp = 1
						if self._data.nmuLoose > 0:
							muPtForIso = max(20., min(self._data.vmuoLoose0_pt, 100.))
							muEtaForIso = min(abs(self._data.vmuoLoose0_eta), 2.3)
							muisoweight = self._muiso_eff.GetBinContent(self._muiso_eff.FindBin(muPtForIso, muEtaForIso))
							muisoweightUp = muisoweight + self._muiso_eff.GetBinError(
								self._muiso_eff.FindBin(muPtForIso, muEtaForIso))
							muisoweightDown = muisoweight - self._muiso_eff.GetBinError(
								self._muiso_eff.FindBin(muPtForIso, muEtaForIso))
							if muisoweight <= 0 or muisoweightDown <= 0 or muisoweightUp <= 0:
								print 'muisoweights are %f, %f, %f, setting all to 1' % (
								muisoweight, muisoweightUp, muisoweightDown)
								muisoweight = 1
								muisoweightDown = 1
								muisoweightUp = 1

						event_weight = pu_weight * k_vjets * k_ttbar * mutrigweight * muidweight * muisoweight
						event_weight_syst = {}
						event_weight_syst["MuTriggerUp"] = pu_weight * k_vjets * k_ttbar * mutrigweightUp * muidweight * muisoweight
						event_weight_syst["MuTriggerDown"] = pu_weight * k_vjets * k_ttbar * mutrigweightDown * muidweight * muisoweight
						event_weight_syst["MuIDUp"] = pu_weight * k_vjets * k_ttbar * mutrigweight * muidweightUp * muisoweight
						event_weight_syst["MuIDDown"] = pu_weight * k_vjets * k_ttbar * mutrigweight * muidweightDown * muisoweight
						event_weight_syst["MuIsoUp"] = pu_weight * k_vjets * k_ttbar * mutrigweight * muidweight * muisoweightUp
						event_weight_syst["MuIsoDown"] = pu_weight * k_vjets * k_ttbar * mutrigweight * muidweight * muisoweightDown
						event_weight_syst["PUUp"] = pu_weight_up * k_vjets * k_ttbar * mutrigweight * muidweight * muisoweight
						event_weight_syst["PUDown"] = pu_weight_down * k_vjets * k_ttbar * mutrigweight * muidweight * muisoweight

				fatjet_pt       = self._data.SelectedJet_pt
				fatjet_eta      = self._data.SelectedJet_eta
				fatjet_msd      = self._data.SelectedJet_msd_puppi
				if self._jet_type == "AK8":
					fatjet_dbtag     = self._data.SelectedJet_doublecsv
				else:
					fatjet_dbtag     = self._data.SelectedJet_doublesub
				fatjet_n2ddt    = self._data.SelectedJet_N2DDT
				fatjet_n2ddt_wp = self._data.SelectedJet_N2DDT_wp
				fatjet_tau21ddt = self._data.SelectedJet_tau21DDT
				fatjet_rho      = self._data.SelectedJet_rho
				fatjet_phi      = self._data.SelectedJet_phi
				fatjet_n2       = self._data.SelectedJet_N2sdb1
				fatjet_nParticles = self._data.SelectedJet_nParticles

				# Inclusive histograms
				self._histograms.GetTH1D("inclusive_pt").Fill(fatjet_pt, event_weight)
				self._histograms.GetTH1D("inclusive_eta").Fill(fatjet_eta, event_weight)
				self._histograms.GetTH1D("inclusive_msd").Fill(fatjet_msd, event_weight)
				self._histograms.GetTH1D("inclusive_n2ddt").Fill(fatjet_n2ddt, event_weight)
				self._histograms.GetTH1D("inclusive_dcsv").Fill(fatjet_dbtag, event_weight)

				# Run selection and fill histograms
				self._event_selectors[selection].process_event(self._data, event_weight)

				# Boxes
				event_boxes = ["all"]
				for wp_pair in self._wps:
					n2ddt_wp = wp_pair[0]
					dbtag_cut = wp_pair[1]
					box_suffix = "_n2wp{}_dbtag{}".format(n2ddt_wp, dbtag_cut)
					if fatjet_n2ddt_wp[n2ddt_wp] <= 0. and fatjet_dbtag >= dbtag_cut:
						event_boxes.append("passn2_passdbtag{}".format(box_suffix))
					if fatjet_n2ddt_wp[n2ddt_wp] <= 0. and fatjet_dbtag < dbtag_cut:
						event_boxes.append("passn2_faildbtag{}".format(box_suffix))
					if fatjet_n2ddt_wp[n2ddt_wp] > 0. and fatjet_dbtag >= dbtag_cut:
						event_boxes.append("failn2_passdbtag{}".format(box_suffix))
					if fatjet_n2ddt_wp[n2ddt_wp] > 0. and fatjet_dbtag < dbtag_cut:
						event_boxes.append("failn2_faildbtag{}".format(box_suffix))

				if self._event_selectors[selection].event_pass():
					for box in event_boxes:
						self._selection_histograms[selection].GetTH2D("{}_pt_vs_msd".format(box)).Fill(fatjet_msd, fatjet_pt, event_weight)
						self._selection_histograms[selection].GetTH2D("{}_pt_vs_msd_unweighted".format(box)).Fill(fatjet_msd, fatjet_pt)
						self._selection_histograms[selection].GetTH1D("{}_nevents".format(box)).Fill(0)
						self._selection_histograms[selection].GetTH1D("{}_nevents_weighted".format(box)).Fill(0, event_weight)
						self._selection_histograms[selection].GetTH1D("{}_pfmet".format(box)).Fill(self._data.pfmet, event_weight)
						self._selection_histograms[selection].GetTH1D("{}_dcsv".format(box)).Fill(fatjet_dbtag, event_weight)
						self._selection_histograms[selection].GetTH1D("{}_n2ddt".format(box)).Fill(fatjet_n2ddt, event_weight)
						self._selection_histograms[selection].GetTH1D("{}_n2".format(box)).Fill(fatjet_n2, event_weight)
						self._selection_histograms[selection].GetTH1D("{}_tau21ddt".format(box)).Fill(fatjet_tau21ddt, event_weight)
						self._selection_histograms[selection].GetTH1D("{}_pt".format(box)).Fill(fatjet_pt, event_weight)
						self._selection_histograms[selection].GetTH1D("{}_msd".format(box)).Fill(fatjet_msd, event_weight)
						self._selection_histograms[selection].GetTH1D("{}_eta".format(box)).Fill(fatjet_eta, event_weight)
						self._selection_histograms[selection].GetTH1D("{}_rho".format(box)).Fill(fatjet_rho, event_weight)
						self._selection_histograms[selection].GetTH2D("{}_met_vs_msd".format(box)).Fill(fatjet_msd, self._data.pfmet, event_weight)

						# Weight systematics: run only if this selection is defined in self._weight_systematics
						if selection in self._weight_systematics:
							for systematic in self._weight_systematics[selection]:
								self._selection_histograms[selection].GetTH2D("{}_{}_pt_vs_msd".format(systematic, box)).Fill(fatjet_msd, fatjet_pt, event_weight_syst[systematic])

					self._selection_histograms[selection].GetTH3D("dbtag_vs_pt_vs_msd").Fill(fatjet_msd, fatjet_pt, fatjet_dbtag)
					self._selection_histograms[selection].GetTH3D("n2ddt_vs_pt_vs_msd").Fill(fatjet_msd, fatjet_pt, fatjet_n2ddt)
					self._selection_histograms[selection].GetTH2D("nparticles_vs_n2").Fill(fatjet_n2, fatjet_nParticles)

	def finish(self):
		if self._output_path == "":
			self._output_path = "/uscms/home/dryu/DAZSLE/data/LimitSetting/InputHistograms_{}.root".format(time.time)
			print "[SignalCutflow::finish] WARNING : Output path was not provided! Saving to {}".format(self._output_path)
		print "[SignalCutflow::finish] INFO : Saving histograms to {}".format(self._output_path)
		f_out = ROOT.TFile(self._output_path, "RECREATE")
		self._histograms.SaveAll(f_out)
		for selection, histogrammer in self._selection_histograms.iteritems():
			histogrammer.SaveAll(f_out)
		for selection, selector in self._event_selectors.iteritems():
			selector.print_cutflow()
			selector.make_cutflow_histograms(f_out)
			selector.save_nm1_histograms(f_out)
		f_out.Close()


if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Produce and plot ieta-iphi histograms to look for buggy events')
	input_group = parser.add_mutually_exclusive_group() 
	input_group.add_argument('--supersamples', type=str, help="List of supersamples to fun, comma separated")
	input_group.add_argument('--files', type=str, help="Input file name(s), comma separated")
	input_group.add_argument('--files_txt', type=str, help="Text files with lists of inputs, comma separated")
	parser.add_argument('--output_folder', type=str, help="Output folder")
	parser.add_argument('--label', type=str, help="Output file label")
	parser.add_argument('--jet_type', type=str, default="AK8", help="AK8 or CA15")
	parser.add_argument('--skim_inputs', action='store_true', help="Run over skim inputs")
	args = parser.parse_args()

	# Make list of input files
	samples = []
	sample_files = {}
	if args.supersamples:
		supersamples = args.supersamples.split(",")
		samples = [] 
		for supersample in supersamples:
			samples.extend(config.samples[supersample])
			for sample in config.samples[supersample]:
				if args.skim_inputs:
					sample_files[sample] = config.skims[sample]
				else:
					sample_files[sample] = config.sklims[sample]
	elif args.files:
		sample_files[args.label] = args.files.split(",")
		samples.append(args.label)
	elif args.files_txt:
		samples.append(args.label)
		sample_files[args.label] = []
		with open(args.files_txt, 'r') as f:
			for line in f:
				sample_files[args.label].append(line.strip())

	for sample in samples:
		print "\n *** Running sample {}".format(sample)
		if "Sbb" in sample or args.skim_inputs or "ZPrime" in sample:
			tree_name = "Events"
		else:
			tree_name = "otree"

		# Sanity check: make sure tree exists in file
		for filename in sample_files[sample]:
			print "[event_selection_histograms] INFO : Checking contents of file {}".format(filename)
			f = ROOT.TFile.Open(filename, "READ")
			t = f.Get(tree_name)
			if not t:
				if tree_name == "otree":
					backup_tree_name = "Events"
				else:
					backup_tree_name = "otree"
				t_backup = f.Get(backup_tree_name)
				if t_backup:
					print "[setup_limits] WARNING : Didn't find tree {} in input file, but did find {}. Changing the tree name, but try to fix this.".format(tree_name, backup_tree_name)
					tree_name = backup_tree_name
				else:
					print "[setup_limits] ERROR : Didn't find tree {} in input file, nor {}. Quitting!".format(tree_name, backup_tree_name)
					sys.exit(1)
			# Check that the "NEvents" histogram is present
			h_NEvents = f.Get("NEvents")
			if not h_NEvents:
				if "data" in sample:
					print "[setup_limits] ERROR : NEvents histogram in not in this file! It is probably corrupt. This is data, so this problem is fatal."
					sys.exit(1)
				else:
					print "[setup_limits] WARNING : NEvents histogram in not in this file! It is probably corrupt. This is MC, so I am skipping the file. But, you probably want to remove from the input list."
					sample_files[sample].remove(filename)
			
		limit_histogrammer = Histograms(sample, tree_name=tree_name, jet_type=args.jet_type)
		output_file_basename ="InputHistograms_{}_{}.root".format(sample, args.jet_type) 
		if args.output_folder:
			limit_histogrammer.set_output_path("{}/{}".format(args.output_folder, output_file_basename))
		else:
			limit_histogrammer.set_output_path("/uscms/home/dryu/DAZSLE/data/LimitSetting/{}".format(output_file_basename))
		for filename in sample_files[sample]:
			print "Input file {}".format(filename)
			limit_histogrammer.add_file(filename)
		#limit_histogrammer.set_jet_type(args.jet_type)
		if "JetHTRun2016" in sample or "SingleMuRun2016" in sample:
			limit_histogrammer.set_data_source("data")
		else:
			limit_histogrammer.set_data_source("simulation")
		if "ps10" in sample:
			limit_histogrammer.set_prescale(10)
		limit_histogrammer.start()
		limit_histogrammer.run()
		limit_histogrammer.finish()

