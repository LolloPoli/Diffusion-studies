#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
//#include <LHAPDF/LHAPDF.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TAxis.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <set>
#include <vector>
#include <TMath.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TPaveStatsEditor.h>
#include "TLorentzVector.h"
#include <cmath>
#include <random>
//#include <TROOT.h>

// per far funzionare il programma
// . /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-centos7-gcc11-opt/setup.sh  intel
// . /cvmfs/sft.cern.ch/lcg/views/LCG_105a/aarch64-el9-gcc13-opt/setup.sh     macOS
// . /cvmfs/sft.cern.ch/lcg/views/LCG_106b/arm64-mac12-clang140-opt/setup.sh  ultima versione
// root
// .L taskB.C
// .L run_all.C
// run()
// NON SO MAGARI SERVONO PER LA NNPDF
// gSystem->AddIncludePath("/cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-centos7-gcc11-opt/include");
//gSystem->AddLinkPath("/cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-centos7-gcc11-opt/lib");

// per non avere le canvas a schermo
gROOT->SetBatch(kTRUE);
gErrorIgnoreLevel = kInfo; // Questo ti mostrerà anche avvisi e warning
gDebug = 0;

std::vector<double> CreateLogBinning(int nbins, double xmin, double xmax) {
    std::vector<double> bin_edges(nbins + 1);
    double logxmin = std::log10(xmin);
    double logxmax = std::log10(xmax);
    double bin_width = (logxmax - logxmin) / nbins;
    for (int i = 0; i <= nbins; ++i) {
        bin_edges[i] = std::pow(10, logxmin + i * bin_width);
    }
    return bin_edges;
}

enum EnPID{
      kelectron,
      kpion,
      kkaon,
      kproton,
      kany
};

/*
double CalcolaMedia(const std::vector<double>& dati) {
    double somma = 0;
    for (double valore : dati) {
        somma += valore;
    }
    return somma / dati.size();
}

double CalcolaErroreStatistico(const std::vector<double>& dati) {
    double media = CalcolaMedia(dati);
    double somma_quad = 0;
    for (double valore : dati) {
        somma_quad += (valore - media) * (valore - media);
    }
    double deviazione_std = sqrt(somma_quad / (dati.size() - 1));
    return deviazione_std / sqrt(dati.size());  // Errore standard della media
}
*/

int checkRecoID(int recpdg){
    switch (std::abs(recpdg)){
        case 0 :
            return kany;
        case 11 : 
            return kelectron;
        case 211 :
            return kpion;
        case 321 :
            return kkaon;
        case 2212 :
            return kproton;
        default :
            Printf ("any particle pdg: %i", recpdg);
            return kany;
    }    
}

double A_siv = 0.02776;
double A_coll = 0.03262;
double PolarizFunction(double Phi_coll, double Phi_siv){
  double funct = 0.5*(1 + A_coll*TMath::Sin(Phi_coll) + A_siv*TMath::Sin(Phi_siv));
  return funct;
}

double RandomProb() {
    static std::random_device rd;                 // Seme casuale (statico, inizializzato una sola volta)
    static std::mt19937 gen(rd());                // Generatore Mersenne Twister
    static std::uniform_real_distribution<> dis(0.0, 1.0); // Distribuzione uniforme

    return dis(gen); // Genera un numero casuale
}

void pino(const char* inputFile1, const char* inputFile2, const char* inputFile3, const char* outputFile){
/*
void pino(TString infile1="pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_4.1469.eicrecon.tree.edm4eic.root",
TString infile2="pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_4.1468.eicrecon.tree.edm4eic.root",
TString infile3="pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_4.1467.eicrecon.tree.edm4eic.root"){
*/  
    // Carica il dizionario per EDM4HEP
    //gSystem->Load("/Users/lorenzopolizzi/Software/EDM4hep/install/lib/libedm4hepDict.dylib");
    // Carica il dizionario per EDM4EIC
    //gSystem->Load("/Users/lorenzopolizzi/Software/edm4eic/install/lib/libedm4eicDict.so");
    gSystem->Load("libedm4hepDict.dylib");
    gSystem->Load("libedm4eicDict.so");    
    // Genera il dizionario mancante
    //gInterpreter->GenerateDictionary("edm4hep::TrackerHitData", "EDM4hepDict.h");
    //Interpreter->GenerateDictionary("edm4hep::ObjectID", "EDM4hepDict.h");

    TFile *ofile = TFile::Open(outputFile, "RECREATE");

    // Set up input file chain
    TChain *mychain = new TChain("events");
    mychain->Add(inputFile1);
    mychain->Add(inputFile2);
    mychain->Add(inputFile3);

    // Initialize reader
    std::vector<std::string> files = {inputFile1, inputFile2, inputFile3};
    for (const auto &file : files) {
        if (mychain->Add(file.c_str()) == 0) {
            std::cerr << "Errore nell'aggiungere il file: " << file << std::endl;
        }
    }
    if (mychain->GetEntries() == 0) {
    std::cerr << "Errore: Nessun evento trovato nella catena!" << std::endl;
    return;
    }
    TTreeReader tree_reader(mychain);
    if (!tree_reader.Next()) {
    std::cerr << "Errore nel lettore di eventi!" << std::endl;
    }
    // Get Particle Information
    TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<float> partMomX(tree_reader, "MCParticles.momentum.x");
    TTreeReaderArray<float> partMomY(tree_reader, "MCParticles.momentum.y");
    TTreeReaderArray<float> partMomZ(tree_reader, "MCParticles.momentum.z");
    TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");
    TTreeReaderArray<int> parentsIndex(tree_reader, "_MCParticles_parents.index");
    TTreeReaderArray<int> daughterIndex(tree_reader, "_MCParticles_daughters.index");
    TTreeReaderArray<unsigned int> par(tree_reader, "MCParticles.parents_end");
    TTreeReaderArray<double> partSpinX(tree_reader, "MCParticles.vertex.x");
    TTreeReaderArray<double> partSpinY(tree_reader, "MCParticles.vertex.y");
    TTreeReaderArray<double> partSpinZ(tree_reader, "MCParticles.vertex.z");
    TTreeReaderArray<double> EndPointX(tree_reader, "MCParticles.endpoint.x");
    TTreeReaderArray<double> EndPointY(tree_reader, "MCParticles.endpoint.y");
    TTreeReaderArray<double> EndPointZ(tree_reader, "MCParticles.endpoint.z");
    TTreeReaderArray<float> dRICHx(tree_reader, "_DRICHAerogelTracks_points.position.x");
    TTreeReaderArray<float> dRICHy(tree_reader, "_DRICHAerogelTracks_points.position.y");
    TTreeReaderArray<float> dRICHz(tree_reader, "_DRICHAerogelTracks_points.position.z");
    TTreeReaderArray<float> dRICH_time(tree_reader, "_DRICHAerogelTracks_points.time");
    TTreeReaderArray<float> dRICH_momX(tree_reader, "_DRICHAerogelTracks_points.momentum.x");
    TTreeReaderArray<float> dRICH_momY(tree_reader, "_DRICHAerogelTracks_points.momentum.y");
    TTreeReaderArray<float> dRICH_momZ(tree_reader, "_DRICHAerogelTracks_points.momentum.z");
    //TTreeReaderArray<float> dRICH_theta(tree_reader, "_DRICHAerogelTracks_points.theta");
    TTreeReaderArray<float> dRICHxGas(tree_reader, "_DRICHGasTracks_points.position.x");
    TTreeReaderArray<float> dRICHyGas(tree_reader, "_DRICHGasTracks_points.position.y");
    TTreeReaderArray<float> dRICHzGas(tree_reader, "_DRICHGasTracks_points.position.z");
    TTreeReaderArray<float> ToF_Rec_x(tree_reader, "TOFEndcapRecHits.position.x");
    TTreeReaderArray<float> ToF_Rec_y(tree_reader, "TOFEndcapRecHits.position.y");
    TTreeReaderArray<float> ToF_Rec_z(tree_reader, "TOFEndcapRecHits.position.z");
    TTreeReaderArray<float> ToF_RecErr_x(tree_reader, "TOFEndcapRecHits.positionError.xx");
    TTreeReaderArray<float> ToF_RecErr_y(tree_reader, "TOFEndcapRecHits.positionError.yy");
    TTreeReaderArray<float> ToF_RecErr_z(tree_reader, "TOFEndcapRecHits.positionError.zz");
    TTreeReaderArray<double> ToFx(tree_reader, "TOFEndcapHits.position.x");
    TTreeReaderArray<double> ToFy(tree_reader, "TOFEndcapHits.position.y");
    TTreeReaderArray<double> ToFz(tree_reader, "TOFEndcapHits.position.z");
    TTreeReaderArray<double> TrackerX(tree_reader, "TrackerEndcapHits.position.x");
    TTreeReaderArray<double> TrackerY(tree_reader, "TrackerEndcapHits.position.y");
    TTreeReaderArray<double> TrackerZ(tree_reader, "TrackerEndcapHits.position.z");
    TTreeReaderArray<float> ForwardMPGD_x(tree_reader, "ForwardMPGDEndcapRecHits.position.x");
    TTreeReaderArray<float> ForwardMPGD_y(tree_reader, "ForwardMPGDEndcapRecHits.position.y");
    TTreeReaderArray<float> ForwardMPGD_z(tree_reader, "ForwardMPGDEndcapRecHits.position.z");
    TTreeReaderArray<float> SiEndcapT_x(tree_reader, "SiEndcapTrackerRecHits.position.x");
    TTreeReaderArray<float> SiEndcapT_y(tree_reader, "SiEndcapTrackerRecHits.position.y");
    TTreeReaderArray<float> SiEndcapT_z(tree_reader, "SiEndcapTrackerRecHits.position.z");

    // Get Reconstructed Track Information
    TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");
    TTreeReaderArray<float> RecEndPointX(tree_reader, "ReconstructedChargedParticles.referencePoint.x");
    TTreeReaderArray<float> RecEndPointY(tree_reader, "ReconstructedChargedParticles.referencePoint.y");
    TTreeReaderArray<float> RecEndPointZ(tree_reader, "ReconstructedChargedParticles.referencePoint.z");
    TTreeReaderArray<int> recPdg(tree_reader, "ReconstructedChargedParticles.PDG");

    // Get Associations Between MCParticles and ReconstructedChargedParticles
    TTreeReaderArray<unsigned int> recoAssoc(tree_reader, "ReconstructedChargedParticleAssociations.recID");
    TTreeReaderArray<unsigned int> simuAssoc(tree_reader, "ReconstructedChargedParticleAssociations.simID");

    
    // GRAFICI DI Q^2
    int nbins = 80;
    int nbon = 60;
    int nben = 80;
    int nbun = 40;
    double xmin_xbj = 1e-2;
    double xmax_xbj = 1;
    double xmin_Q2 = 1;
    double xmax_Q2 = 100.;
    double qm = 1;
    double qM = 1e4;
    double xM = 1;
    std::vector<double> log_bins_Q2 = CreateLogBinning(nbins, xmin_Q2, xmax_Q2);
    std::vector<double> log_bins_xbj = CreateLogBinning(nbins, xmin_xbj, xmax_xbj);
    std::vector<double> log_bins_x = CreateLogBinning(nbon, xmin_xbj, xM);
    std::vector<double> log_bins_Q = CreateLogBinning(nbon, qm, qM);
    
    double xmin_q2 = 1;
    double xmax_q2 = 100;
    double xmin_mom = 1e-1;
    double xmax_mom = 50;
    double xmin_PhT = 1e-2;
    double xmax_PhT = 50;
    double xmin_z = 1e-4;
    double xmax_z = 1;
    std::vector<double> log_bins_Mom = CreateLogBinning(nben, xmin_mom, xmax_mom);
    std::vector<double> log_bins_Q22 = CreateLogBinning(nben, xmin_q2, xmax_q2);
    std::vector<double> log_bins_PhT2 = CreateLogBinning(nben, xmin_PhT, xmax_PhT);
    std::vector<double> log_bins_z2 = CreateLogBinning(nben, xmin_z, xmax_z);
    std::vector<double> log_bins_xbj2 = CreateLogBinning(nben, xmin_xbj, xmax_xbj);
    //TH1D *pion_PhT = new TH1D("Particle_PhT", "Production of Pions in function of P_hT; GeV", nbins,0, 5);
    //TH1D *pion_PT = new TH1D("Particle_PT", "Production of Pions+ in function of P_T (from HP); GeV", nbins, 0, 5);
    TH1D *pion_Q2_dist = new TH1D("pion_Q2_Distribution", "Q^2 distribution | pion| y<0.95 + dRICH acceptance; Q^2 [GeV^2]", nben, log_bins_Q22);
    TH1D *pion_X_dist = new TH1D("pion_Xbj_Distribution", "X_b distribution | pion| y<0.95 + dRICH acceptance; x_b ", nben, log_bins_xbj2);
    TH1D *pion_z_dist = new TH1D("pion_z_Distribution", "X_b distribution | pion| y<0.95 + dRICH acceptance; z ", nben, 0, 1);
    TH1D *pion_PhT_dist = new TH1D("pion_PhT_Distribution", "X_b distribution | pion| y<0.95 + dRICH acceptance; P_hT [GeV] ", nben, 0, 5);
    TH2D *pion_MomVsQ2 = new TH2D("Particle_MomVsQ2", "Mom vs Q^2 | ch. particles | y<0.95; Mom [GeV]; Q^2 [GeV^2]", nben, 2.5, 30, nben, log_bins_Q22.data());
    TH2D *pion_MomVsEta = new TH2D("Particle_MomVsEta", "Mom vs Eta | ch. particles | y<0.95; Eta; Mom [GeV]", nben, 1.5, 3.5, nben, 2.5, 30);
    TH2D *pion_MomVsPhi = new TH2D("Particle_MomVsTheta", "Mom vs Theta (Polar) | ch. particles | y<0.95; Mom [GeV]; Theta [Deg]", nben, 2.5, 30, nben, 3.5, 25);
    TH2D *pion_MomVsTheta = new TH2D("Particle_MomVsPhi", "Mom vs Phi (Azimuth) | ch. particles | y<0.95; Phi [Deg]; Mom [GeV]", nben, -180, 180, nben, 2.5, 30);
    TH2D *pion_MomVsTheta_HP = new TH2D("Particle_MomVsTheta_HP", "Mom vs Theta (Polar, Hadron plane) | ch. particles | y<0.95; Mom [GeV]; Theta [Deg]", nben, 2.5, 30, nben, 3.5, 25);
    TH2D *pion_MomVsPhi_HP = new TH2D("Particle_MomVsPhi_HP", "Mom vs Phi (Azimuth, Hadron plane) | ch. particles | y<0.95; Phi [Deg]; Mom [GeV]", nben, -180, 180, nben, 2.5, 30);
    TH2D *pion_PhTvsQ2 = new TH2D("Particle_PhTVsQ2", "PhT vs Q^2 | ch. particles | y<0.95; PhT [GeV]; Q^2 [GeV^2]", nben, 0, 4, nben, log_bins_Q22.data());
    TH2D *pion_PhTvsMom = new TH2D("Particle_PhTvsMom", "PhT Vs Mom | ch. particles | y<0.95; PhT [GeV]; Mom [GeV]",nben, 0, 4, nben, 2.5, 30);
    TH2D *pion_PhTvsMom2 = new TH2D("Particle_PhTvsMom2", "PhT Vs Mom | ch. particles | y<0.95; PhT [GeV]; Mom [GeV]",40, 0, 1.4, 40, 2.5, 8);
    TH2D *pion_PhTvsZ = new TH2D("Particle_PhTvsZ", "Z vs PhT | ch. particles | y<0.95; z; PhT [GeV]", nben, 0.01, 0.5, nben, 0, 4);
    TH2D *pion_PhTvsEta = new TH2D("Particle_PhTVsEta", "PhT vs Eta | ch. particles | y<0.95; Eta; PhT [GeV]", nben, 1.5, 3.5, nben, 0, 4);
    TH2D *pion_PhTvsPhi = new TH2D("Particle_PhTvsTheta", "Theta (Polar) vs PhT | ch. particles | y<0.95; PhT [GeV]; Theta [Deg]", nben, 0, 4, nben, 3.5, 25);
    TH2D *pion_PhTvsTheta = new TH2D("Particle_PhTvsPhi", "Phi (Azimuth) vs PhT | ch. particles | y<0.95; Phi [Deg]; PhT [GeV]", nben, -180, 180, nben, 0, 4);
    TH2D *pion_PhTvsTheta_HP = new TH2D("Particle_PhTvsTheta_HP", "Theta (Polar, Hadron plane) vs PhT | ch. particles | y<0.95; PhT [GeV]; Theta [Deg]", nben, 0, 4, nben, 3.5, 25);
    TH2D *pion_PhTvsPhi_HP = new TH2D("Particle_PhTvsPhi_HP", "Phi (Azimuth, Hadron plane) vs PhT | ch. particles | y<0.95; Phi [Deg]; PhT [GeV]", nben, -180, 180, nben, 0, 4);
    TH2D *pion_ThetaVsPhi = new TH2D("Particle_PhiVsTheta", "Phi (Azimuth) vs Theta (Polar) | ch. particles | y<0.95; Phi [Deg]; Theta [Deg]", nben, -180, 180, nben, 3.5, 25);
    //TH2D *pion_ThetaVsPhi2 = new TH2D("Particle_PhiVsTheta2", "Phi (Azimuth) vs Theta (Polar) | ch. particles | y<0.95; Phi [Deg]; Theta [Deg]", nben, -180, 180, nben, 0, 20);
    TH2D *pion_Theta_HPVsPhi_HP = new TH2D("Particle_Phi_HPVsTheta_HP", "Phi_HP (Azimuth) vs Theta_HP (Polar)  | ch. particles | y<0.95; Phi [Deg]; Theta [Deg]", nben, -180, 180, nben, 3.5, 25);
    //TH2D *pion_Theta_HPVsPhi_HP2 = new TH2D("Particle_Phi_HPVsTheta_HP2", "Phi_HP (Azimuth) vs Theta_HP (Polar)  | ch. particles | y<0.95; Phi [Deg]; Theta [Deg]", nben, -180, 180, nben, 0, 20);
    TH2D *pion_ThetaVsTheta_HP = new TH2D("Particle_ThetaVsTheta_HP", "Theta (Polar) vs Theta (P, Hadron plane) | ch. particles | y<0.95; Theta [Deg]; Theta [Deg]", nben, 3.5, 25, nben, 3.5, 25);
    TH2D *pion_ThetaVsPhi_HP = new TH2D("Particle_ThetaVsPhi_HP", "Theta (P, Lab) vs Phi (A, Hadron plane) | ch. particles | y<0.95; Phi_HP [Deg]; Theta [Deg]", nben, -180, 180, nben, 3.5, 25);
    //TH2D *pion_ThetaVsPhi_HP2 = new TH2D("Particle_ThetaVsPhi_HP2", "Theta (P, Lab) vs Phi (A, Hadron plane) | ch. particles | y<0.95; Phi_HP [Deg]; Theta [Deg]", nben, -180, 180, nben, 0, 20);
    TH2D *pion_EtaVsPhi_HP = new TH2D("Particle_EtaVsPhi_HP", "Eta vs Phi (A, Hadron plane) | ch. particles | y<0.95; Phi [Deg]; Eta", nben, -180, 180, nben, 1.5, 3.5);
    TH2D *pion_PhiVsPhi_HP = new TH2D("Particle_PhiVsPhi_HP", "Phi (Azimuth) vs Phi (A, Hadron plane) | ch. particles | y<0.95; Phi_lab [Deg]; Phi_h [Deg]", nben, -180, 180, nben, -180, 180);
    TH2D *pion_Phi_HPvsPhi_s = new TH2D("Particle_Phi_HPvsPhi_s", "Phi_h (Hadron plane) vs Phi_s (Spin)  | ch. particles | y<0.95; Phi_h [Deg]; Phi_s [Deg]", nben, -180, 180, nben, -180, 180);
    TH2D *pion_Phi_HPvsPhi_s2 = new TH2D("Particle_Phi_HPvsPhi_s2", "Phi_h + Pi (Hadron plane) vs Phi_s (Spin)  | ch. particles | y<0.95; Phi_h + Pi [Deg]; Phi_s [Deg]", nben, -180, 180, nben, -180, 180);
    TH2D *pion_ZvsQ2 = new TH2D("Particle_ZvsQ2", "Z vs Q^2 | ch. particles | y<0.95; z ; Q^2 [GeV^2]", nben, 0.01, 0.5, nben, log_bins_Q22.data());
    TH2D *pion_ZvsMom = new TH2D("Particle_ZvsMom", "Z vs Mom | ch. particles | y<0.95; z; Mom [GeV]", nben, 0.01, 0.5, nben, 2.5, 30);
    TH2D *pion_ZvsPhi = new TH2D("Particle_ZVsTheta", "Z vs Theta (Polar) | ch. particles | y<0.95; z; Theta [Deg]", nben, 0.01, 0.5, nben, 3.5, 25);
    TH2D *pion_ZvsTheta = new TH2D("Particle_ZvsPhi", "Z vs Phi (Azimuth) | ch. particles | y<0.95; Phi [Deg]; z", nben, -180, 180, nben, 0.01, 0.5);
    TH2D *pion_ZvsTheta_HP = new TH2D("Particle_ZVsTheta_HP", "Z vs Theta (Polar, Hadron plane) | ch. particles | y<0.95; z; Theta [Deg]", nben, 0.01, 0.5, nben, 3.5, 25);
    TH2D *pion_ZvsPhi_HP = new TH2D("Particle_ZvsPhi_HP", "Z vs Phi (Azimuth, Hadron plane) | ch. particles | y<0.95; Phi [Deg]; z", nben, -180, 180, nben, 0.01, 0.5);
    TH2D *pion_P_TvsQ2 = new TH2D("Particle_P_TvsQ2", "P_T vs Q^2 | ch. particles | y<0.95; P_T [GeV]; Q^2 [GeV^2]", nben, 0, 4, nben, log_bins_Q22.data());
    TH2D *pion_P_TvsXb = new TH2D("Particle_P_TvsXb", "Xb vs P_T | ch. particles | y<0.95; xB; PhT [GeV]", nben, log_bins_xbj2.data(), nben, 0, 4);
    TH2D *pion_P_TvsMom = new TH2D("Particle_P_TvsMom", "P_T Vs Mom | ch. particles | y<0.95; P_T [GeV]; Mom [GeV]",nben, 0, 4, nben, 2.5, 30);
    TH2D *pion_P_TvsMom2 = new TH2D("Particle_P_TvsMom2", "P_T Vs Mom | ch. particles | y<0.95; P_T [GeV]; Mom [GeV]",nben, 0, 1.4, nben, 2.5, 8);
    TH2D *pion_P_TvsZ = new TH2D("Particle_P_TvsZ", "Z vs P_T | ch. particles | y<0.95; z; P_T [GeV]", nben, 0.01, 0.5, nben, 0, 4);
    TH2D *pion_P_TvsEta = new TH2D("Particle_P_TVsEta", "P_T vs Eta | ch. particles | y<0.95; Eta; P_T [GeV]", nben, 1.5, 3.5, nben, 0, 4);
    TH2D *pion_P_TvsPhi = new TH2D("Particle_P_TvsTheta", "Theta (Polar) vs P_T | ch. particles | y<0.95; P_T [GeV]; Theta [Deg]", nben, 0, 4, nben, 3.5, 25);
    TH2D *pion_P_TvsTheta = new TH2D("Particle_P_TvsPhi", "Phi (Azimuth) vs P_T | ch. particles | y<0.95; Phi [Deg]; P_T [GeV]", nben, -180, 180, nben, 0, 4);
    TH2D *pion_P_TvsTheta_HP = new TH2D("Particle_P_TvsTheta_HP", "Theta (Polar, Hadron plane) vs P_T | ch. particles | y<0.95; P_T [GeV]; Theta [Deg]", nben, 0, 4, nben, 2.5, 25);
    TH2D *pion_P_TvsPhi_HP = new TH2D("Particle_P_TvsPhi_HP", "Phi (Azimuth, Hadron plane) vs P_T | ch. particles | y<0.95; Phi [Deg]; P_T [GeV]", nben, -180, 180, nben, 0, 4);
    TH2D *pion_Q2vsEta = new TH2D("Particle_Q2vsEta", "Eta vs Q^2 | ch. particles | y<0.95; Eta ; Q^2 [GeV^2]", nben, 1.5, 3.5, nben, log_bins_Q22.data());
    TH2D *pion_XbVsMom = new TH2D("Particle_XbVsMom", "Xb vs Mom | ch. particles | y<0.95; xB; Mom [GeV]", nben, log_bins_xbj2.data(), nben, 2.5, 30);
    TH2D *pion_XbVsPhT = new TH2D("Particle_XbVsPhT", "Xb vs PhT | ch. particles | y<0.95; xB; PhT [GeV]", nben, log_bins_xbj2.data(), nben, 0, 4);
    TH2D *pion_XbVsZ = new TH2D("Particle_XbVsZ", "Xb vs Z | ch. particles | y<0.95; xB; z", nben, log_bins_xbj2.data(), nben, 0.01, 0.5);
    TH2D *pion_XbVsEta = new TH2D("Particle_XbVsEta", "Xb vs Eta | ch. particles | y<0.95; Eta; xB",nben, 1.5, 3.5, nben, log_bins_xbj2.data());
    TH2D *pion_XbVsPhi = new TH2D("Particle_XbVsTheta", "Xb vs Theta (Polar) | ch. particles | y<0.95; xB; Theta [Deg]", nben, log_bins_xbj2.data(), nben, 3.5, 25);
    //TH2D *pion_XbVsPhi2 = new TH2D("Particle_XbVsTheta2", "Xb vs Theta (Polar) | ch. particles | y<0.95; Theta [Deg]; xB", nben, 2.5, 30, nben, log_bins_xbj2.data());
    TH2D *pion_XbVsTheta = new TH2D("Particle_XbvsPhi", "Xb vs Phi (Azimuth) | ch. particles | y<0.95; Phi [Deg]; xB", nben, -180, 180, nben, log_bins_xbj2.data());
    TH2D *pion_XbVsTheta_HP = new TH2D("Particle_XbVsTheta_HP", "Xb vs Theta (Polar, Hadron plane) | ch. particles | y<0.95; xB; Theta [Deg];", nben, log_bins_xbj2.data(), nben, 3.5, 25);
    //TH2D *pion_XbVsTheta_HP2 = new TH2D("Particle_XbVsTheta_HP2", "Xb vs Theta (Polar, Hadron plane) | ch. particles | y<0.95; Theta [Deg]; xB", nben, 2.5, 30, nben, log_bins_xbj2.data());
    TH2D *pion_XbVsPhi_HP = new TH2D("Particle_XbvsPhi_HP", "Xb vs Phi (Azimuth, Hadron plane) | ch. particles | y<0.95; Phi [Deg]; xB", nben, -180, 180, nben, log_bins_xbj2.data());
    TH2D *pion_XbVsQ2 = new TH2D("Particle_XbVsQ2", "Xb vs Q^2 | ch. particles | y<0.95; xB; Q^2 [GeV^2]", nben, log_bins_xbj2.data(), nben, log_bins_Q22.data());
    TH2D *pion_SinPhiVsZ = new TH2D("Particle_SinPhiVsZ", "Sin(Phi_h - Phi_s) vs Z | ch. particles | y<0.95; Sin(); z", nben, -1, 1 , nben, 0.01, 0.5);
    TH2D *pion_SinPhiVsPhT = new TH2D("Particle_SinPhiVsPhT", "Sin(Phi_h - Phi_s) vs PhT | ch. particles | y<0.95; Sin(); PhT [GeV]", nben, -1, 1 , nben, 0, 4);
    TH2D *pion_SinPhiVsMom = new TH2D("Particle_SinPhiVsMom", "Sin(Phi_h - Phi_s) vs Mom | ch. particles | y<0.95; Sin(); Mom [GeV]", nben, -1, 1 , nben, 2.5, 30);
    TH2D *pion_SinPhiVsXb = new TH2D("Particle_SinPhiVsXb", "Sin(Phi_h - Phi_s) vs X_b | ch. particles | y<0.95; Sin(); xB", nben, -1, 1 , nben, log_bins_xbj2.data());
    TH2D *pion_SinPhiVsQ2 = new TH2D("particle_SinPhiVsQ2", "Sin(Phi_h - Phi_s) vs Q2 | ch. particles | y<0.95; Sin(); Q2", nben, -1, 1 , nben, log_bins_Q22.data());    
    TH2D *pion_DeltaPhiVsPhT = new TH2D("Particle_DeltaPhiVsPhT", "Phi_h - Phi_s vs PhT | ch. particles | y<0.95; Phi_h - Phi_s; PhT [GeV]", nben, -180, 180 , nben, 0, 4);
    TH2D *pion_DeltaPhiVsMom = new TH2D("Particle_DeltaPhiVsMom", "Phi_h - Phi_s vs Mom | ch. particles | y<0.95; Phi_h - Phi_s; Mom [GeV]", nben, -180, 180 , nben, 2.5, 30);
    TH2D *pion_DeltaPhiVsXb = new TH2D("Particle_DeltaPhiVsXb", "Phi_h - Phi_s vs Xb | ch. particles | y<0.95; Phi_h - Phi_s; Xb", nben, -180, 180 , nben, log_bins_xbj2.data());
    TH2D *pion_DeltaPhiVsZ = new TH2D("Particle_DeltaPhiVsZ", "Phi_h - Phi_s vs Z | ch. particles | y<0.95; Phi_h - Phi_s; z", nben, -180, 180 , nben, 0.01, 1);
    TH2D *pion_DeltaPhiVsTheta = new TH2D("Particle_DeltaPhiVsTheta", "Phi_h - Phi_s vs Theta | ch. particles | y<0.95; Phi_h - Phi_s; Theta [Deg]", nben, -180, 180 , nben, 3.5, 25);
    TH2D *gamma_ThetaVsPhi2 = new TH2D("gamma_PhiVsTheta2", "Phi (Azimuth) vs Theta (Polar) | gamma | y<0.95; Phi [Deg]; Theta [Deg]", nben, -180, 180, nben, 3.5, 25);

    TH2F* SiTracker_xy = new TH2F("SiTrackerEndcap_xy", "Silicon tracker endcap XY plane; x [mm]; y [mm]", 200, -250, 250, 200, -250, 250);
    TH2D* Tracker_xy = new TH2D("Tracker_xy", "Tracker XY plane distribution ; x [mm]; y [mm]", 200, -250, 250, 200, -250, 250);
    TH2D* Tracker_Forward_xy = new TH2D("Tracker_Forward_xy", "Forward Tracker XY plane distribution ; x [mm]; y [mm]", 200, -250, 250, 200, -250, 250);
    TH2F* MPGD_xy = new TH2F("ForwardMPGD_xy", "Forward MPGD XY plane distribution ; x [mm]; y [mm]", 100, -250, 250, 100, -250, 250);
    TH2D* ToF_xy = new TH2D("ToF_xy", "ToF XY plane distribution ; x [mm]; y [mm]", 100, -250, 250, 100, -250, 250);
    TH2F* ToF_Rec_xy = new TH2F("ToF_Rec_xy", "ToF reconstructed XY plane distribution ; x [mm]; y [mm]", 100, -250, 250, 100, -250, 250);
    TH2F* RealdRICH_xy = new TH2F("RealdRICH_xy", "dRICH XY plane distribution | 1.5<Eta<3.5 , P_h>2.5 GeV ; x [mm]; y [mm]", 100, -250, 250, 100, -250, 250);
    TH2F* RealdRICH_xy_Gas = new TH2F("RealdRICH_xy_Gas", "dRICH XY plane (Gas) distribution | 1.5<Eta<3.5 , P_h>2.5 GeV ; x [mm]; y [mm]", 200, -250, 250, 200, -250, 250);
    TH2D* dRICH_xy = new TH2D("dRICH_xy", "dRICH XY plane distribution | 1.5<Eta<3.5 , P_h>2.5 GeV ; x [mm]; y [mm]", 80, -250, 250, 60, -250, 250);
    TH2F* dRICH_TimeVsMom = new TH2F("dRICH_TimeVsMom", "dRICH Time vs Mom distribution | 1.5<Eta<3.5 , P_h>2.5 GeV ; Mom [GeV]; t [ms]", 150, 0, 20, 150, 1875, 2300);
    //TH2F* dRICH_ThetaVsMom = new TH2F("dRICH_ThetaVsMom", "dRICH Theta vs Mom distribution | 1.5<Eta<3.5 , P_h>2.5 GeV ; Mom [GeV]; Theta [rad]", 150, 0, 20, 150, 0.05, 0.4);

    TH1D *SinPhi_plot = new TH1D("SinPhi", "Sin(Phi_h - Phi_s); Sin()", 50, -1, 1);
    TH1D *Cos2Phi_plot = new TH1D("Cos2Phi_h", "Cos(2*Phi_h); Cos(2Phi_h)", 50, -1, 1);
    TH2D *SinPhivsDelta = new TH2D("SinPhivsDelta", "Sin(Phi_h - Phi_s) vs Phi_h - Phi_s; DeltaPhi; Sin()", 50, -TMath::Pi(), TMath::Pi(), 50, -1, 1);
    TH2D *Cos2PhivsPhi = new TH2D("Cos2Phi_hvsPhi", "Cos(2*Phi_h) vs 2*Phi_h; 2*Phi_h; Cos(2Phi_h)", 50, -2*TMath::Pi(), 2*TMath::Pi(), 50, -1, 1);

    TH1D *P_over_Q = new TH1D("P_t_Over_Q2", "P_T / Q^2; Pt/Q2", 60, 0, 1);
    TH1D *DeltaPhi_plot = new TH1D("DeltaPhi", "Phi_h - Phi_s; Phi_h - Phi_s", 50, -TMath::Pi(), TMath::Pi());
    TH1D *DeltaPhi_plot2 = new TH1D("DeltaPhi2", "Phi_h - Phi_s; Phi_h - Phi_s", 50, -TMath::Pi(), TMath::Pi());
    TH1D *DeltaPhi_plot3 = new TH1D("DeltaPhi_Sin", "Sivers | Phi_h - Phi_s (w: Sin(delta)); Phi_h - Phi_s", 40, -TMath::Pi(), TMath::Pi());
    TH1D *DeltaPhi_sinArea = new TH1D("DeltaPhi_Sin_Area", "Sivers | Phi_h - Phi_s (area); Phi_h - Phi_s", 40, -TMath::Pi(), TMath::Pi());
    TH1D *DeltaPhi_xLow  = new TH1D("DeltaPhi_xLow", "Sivers | Phi_h - Phi_s (w: Sin(delta)) | 0.01 < xb < 0.1; Phi_h - Phi_s", 40, -TMath::Pi(), TMath::Pi());
    TH1D *DeltaPhi_xHigh  = new TH1D("DeltaPhi_xHigh", "Sivers | Phi_h - Phi_s (w: Sin(delta)) | xb > 0.1; Phi_h - Phi_s", 40, -TMath::Pi(), TMath::Pi());
    TH1D *DeltaPhi_Polarized = new TH1D("DeltaPhi_Polarized", "Sivers | Phi_h - Phi_s (w: Sin(delta)); Phi_h - Phi_s", 40, -TMath::Pi(), TMath::Pi());
    TH1D *DeltaPhi_Polarized_Area = new TH1D("DeltaPhi_Polarized_Area", "Sivers | Phi_h - Phi_s (w: Sin(delta)); Phi_h - Phi_s", 40, -TMath::Pi(), TMath::Pi());
    TH1D *DeltaPhi_Pol_Collins = new TH1D("DeltaPhi_Pol_Collins", "Collins | Phi_h - Phi_s (w: Sin(delta)); Phi_h - Phi_s", 40, -TMath::Pi(), TMath::Pi());
    TH1D *DeltaPhi_Pol_CollinsArea = new TH1D("DeltaPhi_Pol_CollinsArea", "Collins | Phi_h - Phi_s (w: Sin(delta)); Phi_h - Phi_s", 40, -TMath::Pi(), TMath::Pi());
    TH1D *DeltaPhi_Collins = new TH1D("DeltaPhi_Collins", "Collins | Phi_h + Phi_s (w: Sin(delta)); Phi_h + Phi_s", 40, -TMath::Pi(), TMath::Pi());
    TH1D *DeltaPhi_CollinsArea = new TH1D("DeltaPhi_Collins_Area", "Collins | Phi_h + Phi_s (area); Phi_h + Phi_s", 40, -TMath::Pi(), TMath::Pi());

    TH1D *Beam_Angle = new TH1D("Beam_Angle", "Angle between the beam; Deg", 60, 178, 179);
    TH1D *Phi_h = new TH1D("Phi_h", "Phi_h; Phi_h [Deg]", 60, -180, 180);
    TH1D *Phi_h2 = new TH1D("Phi_h2", "Phi_h Abs; Phi_h [Deg]", 60, 0, 180);
    TH1D *Phi_ss = new TH1D("Phi_s", "Phi_s; Phi_s [Deg]", 60, -180, 180);
    TH1D *Phi_ss2 = new TH1D("Phi_s2", "Phi_s Abs; Phi_s [Deg]", 60, 0, 180);
    TH1D *z_plot = new TH1D("z_plot", "Z; z", nben, 0.01, 0.5);

    // RICOSTRUZIONI ___________________________________________________________________________________________________________________________________________________________________________
    // ______________________________________________________________________________________________________________________________________________

    // CERCA DEI PROTONI LANCIATI
    double xmin_m = 1e-1;
    double xmax_m = 300;
    double xmax_m2 = 20;
    std::vector<double> log_bins_momA = CreateLogBinning(nbins, xmin_m, xmax_m);
    std::vector<double> log_bins_momB = CreateLogBinning(nbins, xmin_m, xmax_m2);
    TH1D *beamproton_mom = new TH1D("BeamProton_Mom", "Proton beam; GeV",  80, log_bins_momA.data());
    TH1D *beamproton_theta = new TH1D("BeamProton_Phi", "Proton beam angle", 80, 0, 360);
    TH1D *beamproton_Phi = new TH1D("BeamProton_Theta", "Proton beam angle", 80, 0, 360);
    TH1D *beamelectron_mom = new TH1D("BeamElectron_Mom", "Electron beam; GeV",  80, log_bins_momB.data());
    TH1D *beamelectron_theta = new TH1D("BeamElectron_Phi", "Electron beam angle", 80, 0 , 360);
    TH1D *beamelectron_Phi = new TH1D("BeamElectron_Theta", "Electron beam angle", 80, 0 , 360);
    TGraph2D *graphProton = new TGraph2D();      // Grafico per i protoni
    TGraph2D *graphElectron = new TGraph2D();    // Grafico per gli elettroni
    int protonIndex = 0, electronIndex = 0; 
    TH2D *h2 = new TH2D("h2", "Proiezione sferica; Theta (deg); Phi (deg); Intensità", 
                        20, 3.5, 25, 20, -180, 180);
    TH2D *h3 = new TH2D("h3", "Proiezione sferica; Theta (deg); Phi (deg); Intensità", 
                        20, 3.5, 25, 20, -180, 180);
    //PROVA DI TH3D
    /*
    int nban = 10;
    std::vector<double> log_bins_Mom3 = CreateLogBinning(nban, xmin_mom, xmax_mom);
    std::vector<double> log_bins_PhT3 = CreateLogBinning(nban, xmin_PhT, xmax_PhT);
    std::vector<double> log_bins_z3 = CreateLogBinning(nban, xmin_z, xmax_z);
    std::vector<double> log_bins_xbj3 = CreateLogBinning(nban, xmin_xbj, xmax_xbj);
    TH3D *pion_MomVsPhiVsTheta = new TH3D("Particle_MomVsPhiVsTheta", 
    "Mom vs Theta (Polar) vs Phi (Azimuth) | ch. particles | y<0.95; Phi [Deg]; Mom [GeV]; Theta [Deg]", 
    nban, -180, 180, nban, 2.5, 30, nban, 3.5, 25);
    pion_MomVsPhiVsTheta->GetYaxis()->Set(nban, log_bins_Mom3.data());
    TH3D *pion_PhTvsZvsXb = new TH3D("Particle_PhTvsZvsXb",
    "PhT vs Z vs xB | ch. particles | y<0.95; xB; z; PhT [GeV]",
    nban, 1e-2, 1, nban, 0.01, 0.5, nban, 0, 4);
    pion_PhTvsZvsXb->GetXaxis()->Set(nban, log_bins_xbj3.data());
    pion_PhTvsZvsXb->GetYaxis()->Set(nban, log_bins_z3.data());
    pion_PhTvsZvsXb->GetZaxis()->Set(nban, log_bins_PhT3.data());
    TH3D *z_axisSP = new TH3D("z_axisSP", "Scattered photon polar angle; z; x; y",  
    30, -1, 1, 30, -1, 1, 30, -1, 1);
    TH3D *axisHP = new TH3D("Particle_axisSP", "Particle in the scattered plane system; z; x; y",  
    30, -1, 1, 30, -1, 1, 30, -1, 1);
    TH3D *Vertex_Proton = new TH3D("Vertex_proton", "Vertex proton; z; x; y",  
    30, -1, 1, 30, -1, 1, 30, -1, 1);
    TH3D *Vertex_electron = new TH3D("Vertex_electron", "Vertex electron; z; x; y",  
    30, -1, 1, 30, -1, 1, 30, -1, 1);
    TH3D *Vertex_gamma = new TH3D("Vertex_gamma", "Vertex gamma; z; x; y",  
    30, -1, 1, 30, -1, 1, 30, -1, 1);
    TH3D *Vertex_gamma2 = new TH3D("Vertex_gamma2", "Vertex gamma; z; x; y",  
    30, 0.9, 1, 30, -0.15, 0.15, 30, -0.2, 0.2);
    TH3D *Vertex_pion_Unit = new TH3D("Vertex_pion_Unit", "Vertex pion; z; x; y",  
    30, -1, 1, 30, -1, 1, 30, -1, 1;
    TH3D *Vertex_pion = new TH3D("Vertex_pion", "Vertex pion; Pz [GeV]; Px [GeV]; Py [GeV]",  
    20, -2, 20, 20, -0.6, 0.6, 20, -0.4, 0.4);
    */

    /*
    // VEDIAMO SE RIUSCIAMO A FARE LE PDF
    LHAPDF::PDFSet pdfSet("NNPDF31_nnlo_as_0118");
    auto pdf = pdfSet.mkPDF(0);  // '0' è l'indice della replica centrale
    TH2F* hPDF = new TH2F("hPDF", "PDF u(x, Q^{2})", nben, log_bins_xbj2.data(), nben, log_bins_Q22.data());
    hPDF->GetXaxis()->SetTitle("x");
    hPDF->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
    hPDF->GetZaxis()->SetTitle("x f_{u}(x, Q^{2})");
    int pion_count = 0;
    TGraphErrors *uPDF = new TGraphErrors();
    TGraphErrors *dPDF = new TGraphErrors();
    TGraphErrors *glPDF = new TGraphErrors();
    std::vector<double> val_uPDF;
    std::vector<double> val_dPDF;
    std::vector<double> val_glPDF;
    std::vector<double> val_x;
    */
    // FOTONEEE
    TH1D *GammaPolarDeg = new TH1D("GammaPolarDeg", "Scattered photon polar angle; Theta [Deg]",  nben, 3.5, 25);
    TH1D *GammaAzimuthDeg = new TH1D("GammaAzimuthDeg", "Scattered photon azimuthal angle; Phi [Deg]",  nben, -180, 180);


    // ALCUNI VETTORI UTILI
    std::vector<TLorentzVector> scatElectron;
    std::vector<TVector3> recScatElectron;
    std::vector<TVector3> GammaVector;
    std::vector<TVector3> BeamElectronVector;
    TVector3 GammaVectorSingle;
    TVector3 BeamElectronVectorSingle;
    TVector3 GammaVectorSingle_Rec;
    TVector3 BeamElectronVectorSingle_Rec;
    TVector3 BeamPr;
    TVector3 ProtonSpin(0, 1, 0);
    TVector3 ProtonSpinDown(0, -1, 0);
    TVector3 z_lab(0,0,1);
    //std::vector<TVector3> z_axis_SP 
    //std::vector<TVector3> ipsilon;
    //std::vector<TVector3> y_axis_SP;
    //std::vector<TVector3> x_axis_SP;
    std::vector<double> scatPhi;
    std::vector<double> PolarGammaRad;
    std::vector<double> PolarGammaDeg;
    std::vector<double> AzimuthGammaRad;
    std::vector<double> AzimuthGammaDeg;
    std::vector<TVector3> elMom_pion;
    TVector3 ElBeam(0.,0.,-18.); 
    TLorentzVector ElectronBeam;
    TLorentzVector ElectronScattered;
    TLorentzVector ProtonBeam;
    TLorentzVector VirtualPhoton;
    std::vector<TLorentzVector> q;
    double currentPhi = 0;
    double currentMom = 0;
    std::vector<double> scatElPhipion;
    // per la costruzione del pione
    TVector3 currentQ2pion;
    std::vector<TVector3> scatElq_pion;
    // del kaone
    TVector3 currentQ2kaon;
    std::vector<TVector3> scatElq_kaon;
    // del protone
    TVector3 currentQ2proton;
    std::vector<TVector3> scatElq_proton;
    double count = 0;
    double countEta = 0;
    int count_el = 0;
    std::set<int> uniqueStatuses;
    std::set<int> uniqueParentsIndex;
    float dRICH_x, dRICH_y;
    // COUNT PARTICLE SPIN
    double spinUp, spinDown;
    float dRICH_xGas, dRICH_yGas;
    double Phi_s_Down, Phi_s_Up, Delta_Down, Delta_Up, SinPhiDown, SinPhiUp;
    double Delta_Down_Coll, Delta_Up_Coll, SinPhiDown_Coll, SinPhiUp_Coll;
    float ToF_Rx, ToF_Ry, ToF_ErrX, ToF_ErrY;
    double ToF_x, ToF_y;
    double Tracker_x, Tracker_y, Tracker_z;
    float mpgd_x, mpgd_y;
    float SiEnd_x, SiEnd_y, SiEnd_z;
    float drich_time, drich_momX, drich_momY, drich_momZ, dRICH_momentum, drich_theta;

    while(tree_reader.Next()) { // Loop over events
      for(unsigned int l=0; l<dRICHx.GetSize(); l++){
        dRICH_x = dRICHx[l];
        dRICH_y = dRICHy[l];
        drich_time = dRICH_time[l];
        TVector3 dRICH_Mom(dRICH_momX[l],dRICH_momY[l],dRICH_momZ[l]);
        dRICH_momentum = dRICH_Mom.Mag();
        //drich_theta = dRICH_theta[l];
        RealdRICH_xy->Fill(dRICH_x, dRICH_y);
        dRICH_TimeVsMom->Fill(dRICH_momentum, drich_time);
        //dRICH_ThetaVsMom->Fill(dRICH_momentum, drich_theta);
      }
      for(unsigned int a=0; a<dRICHxGas.GetSize(); a++){
        dRICH_xGas = dRICHxGas[a];
        dRICH_yGas = dRICHyGas[a];
        RealdRICH_xy_Gas->Fill(dRICH_xGas, dRICH_yGas);
      }
      for(unsigned int b=0; b<ToF_Rec_x.GetSize(); b++){
        ToF_Rx = ToF_Rec_x[b];
        ToF_Ry = ToF_Rec_y[b];
        ToF_ErrX = ToF_RecErr_x[b];
        ToF_ErrY = ToF_RecErr_y[b];
        double err = sqrt(pow(ToF_ErrX, 2) + pow(ToF_ErrY, 2));
        ToF_Rec_xy->Fill(ToF_Rx, ToF_Ry, err);
      }
      for(unsigned int c=0; c<ToFx.GetSize(); c++){
        ToF_x = ToFx[c];
        ToF_y = ToFy[c];
        ToF_xy->Fill(ToF_x, ToF_y);
      }
      for(unsigned int d=0; d<TrackerX.GetSize(); d++){
        Tracker_x = TrackerX[d];
        Tracker_y = TrackerY[d];
        Tracker_z = TrackerZ[d];
        Tracker_xy->Fill(Tracker_x, Tracker_y);
        if(Tracker_z > 0){
          Tracker_Forward_xy->Fill(Tracker_x, Tracker_y);
        }
      }
      for(unsigned int e=0; e<ForwardMPGD_x.GetSize(); e++){
        mpgd_x = ForwardMPGD_x[e];
        mpgd_y = ForwardMPGD_y[e];
        MPGD_xy->Fill(mpgd_x, mpgd_y);
      }
      for(unsigned int f=0; f<SiEndcapT_x.GetSize(); f++){
        SiEnd_x = SiEndcapT_x[f];
        SiEnd_y = SiEndcapT_y[f];
        SiEnd_z = SiEndcapT_z[f];
        if (SiEnd_z > 0){
          SiTracker_xy->Fill(SiEnd_x, SiEnd_y);
        }
      }

      for(unsigned int i=0; i<partGenStat.GetSize(); i++) // Loop over thrown particles
        {
          count += 1;
          TVector3 part(partMomX[i],partMomY[i],partMomZ[i]);
          double partMom = part.Mag();
          double partEta = part.PseudoRapidity();
          double phis = part.Theta();
          double y_DA= ((TMath::Tan(phis*0.5)) / (TMath::Tan(phis*0.5) + TMath::Tan(currentPhi*0.5)));
          if(partEta >= 1.5 && partEta <= 3.5){
            if(partMom >= 2.5){
                if(y_DA <= 0.95){
                countEta += 1;
                }
            }
          }

          int pdg = (std::abs(partPdg[i]));
          // status = 4 is the beam (ref 1767) in HepMC
          /*
          if(pdg == 11){
           uniqueStatuses.insert(parentsIndex[i]);
          }
          */
          // status = 21 incoming particles of the hardest subprocess 
        
          if(partGenStat[i] == 4){
            if(pdg == 2212){
                TVector3 BeamProton(partMomX[i],partMomY[i],partMomZ[i]);
                TVector3 VertexPr(partSpinX[i],partSpinY[i],partSpinZ[i]);
                BeamPr.SetXYZ(partMomX[i],partMomY[i],partMomZ[i]);
                TVector3 SponPr = BeamProton.Unit();
                double prot_vX = SponPr.X();
                double prot_vY = SponPr.Y();
                double prot_vZ = SponPr.Z();
                //Vertex_Proton->Fill(prot_vZ, prot_vX, prot_vY);
                double Pmom = BeamProton.Mag();
                double Peng = sqrt(Pmom*Pmom + 0.939*0.939);
                double Beam_pThetaRad = BeamProton.Phi();
                double Beam_pThetaDeg = Beam_pThetaRad * (180.0/TMath::Pi());
                double Beam_pPhiRad = BeamProton.Theta();
                double Beam_pPhiDeg = Beam_pPhiRad * (180.0/TMath::Pi());
                double Px = std::abs(1) * TMath::Sin(Beam_pPhiRad) * TMath::Cos(Beam_pThetaRad);
                double Py = std::abs(1) * TMath::Sin(Beam_pPhiRad) * TMath::Sin(Beam_pThetaRad);
                double Pz = std::abs(1) * TMath::Cos(Beam_pPhiRad);

                graphProton->SetPoint(protonIndex++, Px, Py, Pz);
                beamproton_mom->Fill(Peng);
                beamproton_theta->Fill(Beam_pThetaDeg);
                beamproton_Phi->Fill(Beam_pPhiDeg);
                h2->Fill(Beam_pPhiDeg, Beam_pThetaDeg, 1);
                ProtonBeam.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i],Pmom);
              }
              if(pdg == 11){
                TVector3 BeamElectron(partMomX[i],partMomY[i],partMomZ[i]);
                TVector3 VertexE(partSpinX[i],partSpinY[i],partSpinZ[i]);
                TVector3 SponE = BeamElectron.Unit();
                double el_vX = SponE.X();
                double el_vY = SponE.Y();
                double el_vZ = SponE.Z();
                //Vertex_electron->Fill(el_vZ, el_vX, el_vY);
                double Emom = BeamElectron.Mag();
                double Beam_EThetaRad = BeamElectron.Phi();
                double Beam_EThetaDeg = Beam_EThetaRad * (180.0/TMath::Pi());
                double Beam_EPhiRad = BeamElectron.Theta();
                double Beam_EPhiDeg = Beam_EPhiRad * (180.0/TMath::Pi());
                double Ex = std::abs(1) * TMath::Sin(Beam_EPhiRad) * TMath::Cos(Beam_EThetaRad);
                double Ey = std::abs(1) * TMath::Sin(Beam_EPhiRad) * TMath::Sin(Beam_EThetaRad);
                double Ez = std::abs(1) * TMath::Cos(Beam_EPhiRad);
                double BeamAngle = BeamElectron.Angle(BeamPr);
                double BeamAngle_Deg = BeamAngle * (180 / TMath::Pi());
                Beam_Angle->Fill(BeamAngle_Deg);
                ElectronBeam.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i],Emom);

                graphElectron->SetPoint(electronIndex++, Ex, Ey, Ez);
                beamelectron_mom->Fill(Emom);
                beamelectron_theta->Fill(Beam_EThetaDeg);
                beamelectron_Phi->Fill(Beam_EPhiDeg);
                h3->Fill(Beam_EPhiDeg, Beam_EThetaDeg, 1);
              }
            }
          if(partGenStat[i] <= 1)
            {
              if(pdg == 11)
                {
                  TVector3 ElMom(partMomX[i],partMomY[i],partMomZ[i]);
                  double mom = ElMom.Mag();
                  double etta = ElMom.PseudoRapidity();
                  if(std::abs(etta) <= 3.5){
                    if(parentsIndex[i]<=5){
                      TLorentzVector tlv(mom, ElMom.X(), ElMom.Y(), ElMom.Z());
                      double angleR = ElMom.Theta();
                      double angle = angleR * (180.0 / TMath::Pi());
                      currentPhi = ElMom.Theta();
                      currentMom = ElMom.Mag();
                      currentQ2pion.SetXYZ(partMomX[i],partMomY[i],partMomZ[i]);
                      currentQ2kaon.SetXYZ(partMomX[i],partMomY[i],partMomZ[i]);
                      currentQ2proton.SetXYZ(partMomX[i],partMomY[i],partMomZ[i]);
                      //scatElectron.push_back(ElMom);  // to use it outside the cycle
                      scatElectron.push_back(tlv);
                      scatPhi.push_back(angle);
                      //BeamElectronVector.push_back(ElMom);
                      BeamElectronVectorSingle.SetXYZ(partMomX[i],partMomY[i],partMomZ[i]);
                      ElectronScattered.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i], currentMom);

                        for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
                        {
                          if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
                            {
                              TVector3 recElmom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); 
                              double momE = recElmom.Mag();
                              //recScatEl->Fill(momE);
                              recScatElectron.push_back(recElmom);
                              BeamElectronVectorSingle_Rec.SetXYZ(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); 
                            }
                        }
                    }
                  }
                }
                if(pdg == 22){ // stiamo cercando il fotone scatterato
                  if(parentsIndex[i]<=6){
                    TVector3 ScatPhoton(partMomX[i],partMomY[i],partMomZ[i]);
                    TVector3 VertexG(partSpinX[i],partSpinY[i],partSpinZ[i]);
                    double GammaMom = ScatPhoton.Mag();
                    TVector3 SponG = ScatPhoton.Unit();
                    double gamma_vX = SponG.X();
                    double gamma_vY = SponG.Y();
                    double gamma_vZ = SponG.Z();
                    //Vertex_gamma->Fill(gamma_vZ, gamma_vX, gamma_vY);
                    //Vertex_gamma2->Fill(gamma_vZ, gamma_vX, gamma_vY);
                    double Polar_gammaRad = ScatPhoton.Theta();
                    double Polar_gammaDeg = Polar_gammaRad * (180.0/TMath::Pi());
                    double Azimuth_gammaRad = ScatPhoton.Phi();
                    double Azimuth_gammaDeg = Azimuth_gammaRad * (180.0/TMath::Pi());
                    PolarGammaDeg.push_back(Polar_gammaDeg);
                    PolarGammaRad.push_back(Polar_gammaRad);
                    AzimuthGammaDeg.push_back(Azimuth_gammaDeg);
                    AzimuthGammaRad.push_back(Azimuth_gammaRad);
                    GammaPolarDeg->Fill(Polar_gammaDeg);
                    GammaAzimuthDeg->Fill(Azimuth_gammaDeg);
                    gamma_ThetaVsPhi2->Fill(Azimuth_gammaDeg, Polar_gammaDeg);
                    //GammaVector.push_back(ScatPhoton);
                    GammaVectorSingle.SetXYZ(partMomX[i],partMomY[i],partMomZ[i]);
                    VirtualPhoton.SetPxPyPzE(partMomX[i],partMomY[i],partMomZ[i],GammaMom);
                    for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
                        {
                          if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
                            {
                              TVector3 ScatPhoton_Rec(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); 
                              TVector3 SponG_Rec = ScatPhoton_Rec.Unit();
                              double gamma_vX_Rec = SponG_Rec.X();
                              double gamma_vY_Rec = SponG_Rec.Y();
                              double gamma_vZ_Rec = SponG_Rec.Z();
                              GammaVectorSingle_Rec.SetXYZ(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
                              
                            }
                        }
                  }
                }
            }
          // IDK IF THIS WILL WORK HERE...
          for (const auto& vec : scatElectron) {
            q.push_back(ElectronBeam - vec);
          }
          if(partGenStat[i] == 1) // Select stable thrown particles
            {  
              if(pdg == 211 || pdg == 321 || pdg == 2212)
                { 
                  if(parentsIndex[i]<=10){
                      TVector3 truePionMom(partMomX[i],partMomY[i],partMomZ[i]);
                      TVector3 VertexPi(partSpinX[i],partSpinY[i],partSpinZ[i]);
                      TVector3 EndPointVector(EndPointX[i],EndPointY[i],EndPointZ[i]);
                      //TVector3 dRICH_position(dRICHx[i],dRICHy[i],dRICHz[i]);
                      TVector3 AxPi = truePionMom.Unit();
                      double pion_vX = truePionMom.X();
                      double pion_vY = truePionMom.Y();
                      double pion_vZ = truePionMom.Z();
                      double pion_vX_Unit = AxPi.X();
                      double pion_vY_Unit = AxPi.Y();
                      double pion_vZ_Unit = AxPi.Z();
                      // PIANO DEL dRICH
                      double pion_EpZ = EndPointVector.Z();
                      double t = (1980 - pion_EpZ) / pion_vZ_Unit;
                      double x_drich = EndPointVector.X() + t*pion_vX_Unit;
                      double y_drich = EndPointVector.Y() + t*pion_vY_Unit;
                      //Vertex_pion->Fill(0.1*pion_vZ, 0.1*pion_vX, 0.1*pion_vY);
                      //Vertex_pion_Unit->Fill(pion_vZ_Unit, pion_vX_Unit, pion_vY_Unit);
                      double mom_pion = truePionMom.Mag();
                      double E_pion = sqrt(mom_pion*mom_pion + 0.139*0.139);
                      TLorentzVector pion(partMomX[i],partMomY[i],partMomZ[i], mom_pion);
                      double pionEta = truePionMom.PseudoRapidity();

                      if(pionEta >= 1.5 && pionEta <= 3.5){
                        if(mom_pion >= 2.5){
                          double pionPhiRad = truePionMom.Theta(); // MOLTO CONFUSO LO SO, PHI (MIO) ANGOLO POLARE
                          double pionThetaRad = truePionMom.Phi();
                          /*
                          double pionThetaRad = truePionMom.Phi() + TMath::Pi(); // THETA ANGOLO AZIMUTALE 
                          if(pionThetaRad >= TMath::Pi()){
                            pionThetaRad -= 2*TMath::Pi();
                          }*/
                          double pionPhi = pionPhiRad * (180.0 / TMath::Pi());
                          double pionTheta = pionThetaRad * (180.0 / TMath::Pi());
                          TLorentzVector scatElpion(currentMom, currentQ2pion.X(), currentQ2pion.Y(), currentQ2pion.Z());
                          //TLorentzVector photon_pion = ElectronBeam - scatElpion; 
                          TLorentzVector photon_pion = ElectronBeam - ElectronScattered;
                          scatElPhipion.push_back(currentPhi); // that generate an array with the angle of the scat. el with the # of pions
                          scatElq_pion.push_back(currentQ2pion);
                          double y_DA_pion = ((TMath::Tan(pionPhiRad*0.5)) / (TMath::Tan(pionPhiRad*0.5) + TMath::Tan(currentPhi*0.5)));
                          if(y_DA_pion <= 0.95){  
                              //TVector3 CurrentGammaVector = GammaVector[i];
                              TVector3 CurrentGammaVector = GammaVectorSingle;
                              TVector3 z_axis_SP = CurrentGammaVector.Unit(); // voglio che l'asse z sia lungo la direzione di gamma
                              double z_axis_z = z_axis_SP.Z();
                              double z_axis_y = z_axis_SP.Y();
                              double z_axis_x = z_axis_SP.X();
                              // per calcolare y mi serve il prodotto vettoriale tra il leptone entrante e gamma
                              TVector3 CurrentBeamElectronVector = BeamElectronVectorSingle;
                              TVector3 ipsilon = CurrentGammaVector.Cross(CurrentBeamElectronVector);
                              TVector3 y_axis_SP = ipsilon.Unit();
                              TVector3 x_axis_SP = y_axis_SP.Cross(z_axis_SP);
                              // BENE ADESSO ABBIAMO I NOSTRI ASSI DEL PIANO DI SCATTERING
                              double P_T = truePionMom.Perp(z_axis_SP); // questo comando mi da il vettore ortogonale a z, quindi il momento trasverso
                              double Theta_HP_Rad = truePionMom.Angle(z_axis_SP); // dovrebbe fornire l'angolo polare del hadron plane
                              double Theta_HP_Deg = Theta_HP_Rad * (180.0/TMath::Pi());
                              // però vogliamo anche il vettore del momento trasverso
                              double PdotZ = truePionMom.Dot(z_axis_SP); // momento proiettato sull'asse z
                              double PdotY = truePionMom.Dot(y_axis_SP);
                              double PdotX = truePionMom.Dot(x_axis_SP);
                              TVector3 axis_HP(PdotX, PdotY, PdotZ);
                              TVector3 axis_HP_Unit = axis_HP.Unit();
                              double z = axis_HP_Unit.Z();
                              double y = axis_HP_Unit.Y();
                              double x = axis_HP_Unit.X();
                              TVector3 Momentum_Z = PdotZ * z_axis_SP;
                              TVector3 P_T_Vector = truePionMom - Momentum_Z; // ora lo possiamo usare per trovare l'angolo azimutale
                              TVector3 P_T_axis = P_T_Vector.Unit();
                              double P_T_Vector_x = P_T_Vector * x_axis_SP;
                              double P_T_Vector_y = P_T_Vector * y_axis_SP; // calcolo il momento trasverso sugli assi x e y
                              double Phi_HP_Rad = std::atan2(P_T_Vector_y, P_T_Vector_x);
                              double Phi_HP_Deg = Phi_HP_Rad * (180.0/TMath::Pi()); 
                              double phi_hh_Rad = P_T_Vector.Angle(x_axis_SP);
                              double phi_hh_Deg = phi_hh_Rad * (180/TMath::Pi());
                              double Phi_acca = Phi_HP_Rad + TMath::Pi();
                              if(Phi_acca >= TMath::Pi()){
                                Phi_acca -= 2*TMath::Pi();
                              }
                              double Phi_acca_Deg = Phi_acca * (180.0/TMath::Pi()); 
                              // MISURA DI PHI_S
                              //double Phi_s = x_axis_SP.Angle(ProtonSpin);
                              TVector3 termine1_1 = x_axis_SP.Cross(ProtonSpin);
                              double termine1 = termine1_1 * z_axis_SP;
                              double termine2 = x_axis_SP.Dot(ProtonSpin);
                              double Phi_s = std::atan2(termine1, termine2);
                              double Phi_s_Deg = Phi_s * (180/TMath::Pi());
                              double phi_s22 = x_axis_SP.Angle(ProtonSpin);
                              double phi_s22_Deg = phi_s22 * (180/TMath::Pi());
                              // CALCOLO DIPENDENZE SINUSOIDALI
                              double DeltaPhi = Phi_HP_Rad - Phi_s;
                              double DeltaPhi2 = phi_hh_Rad - phi_s22;
                              if(DeltaPhi <-TMath::Pi()){
                                DeltaPhi += 2*TMath::Pi();
                              }
                              else if(DeltaPhi > TMath::Pi()){
                                DeltaPhi -= 2*TMath::Pi();
                              }
                              double DeltaPhi_Rad = DeltaPhi * (180.0/TMath::Pi());
                              double DeltaPhi_Rad2 = DeltaPhi2 * (180.0/TMath::Pi());
                              double SinPhi = TMath::Sin(DeltaPhi);
                              double Cos2Phi_h = TMath::Cos(2*Phi_HP_Rad);
                              double Q2_DA_pion = ((4*18*18*(1-y_DA_pion)) / (TMath::Tan(currentPhi*0.5)*TMath::Tan(currentPhi*0.5)));
                              //double Q2_DA_pion = (scp*scp*TMath::Sin(currentPhi)*TMath::Sin(currentPhi))/(1-y_DA_pion);
                              double xbj_DA_pion = (Q2_DA_pion) / (4*18*275*y_DA_pion);
                              double z_DA_pion = (ProtonBeam * pion) / (ProtonBeam * photon_pion);
                              TVector3 PhT_pion_vec = truePionMom - ((truePionMom * currentQ2pion) / currentQ2pion.Mag())*currentQ2pion.Unit();
                              double PhT_pion = PhT_pion_vec.Mag();
                              TLorentzVector sq = ElectronBeam + ProtonBeam;
                              double Rapporto = P_T / Q2_DA_pion;
                              // dRICH position
                              //dRICH_x = dRICH_position.X();
                              //dRICH_y = dRICH_position.Y();
                              // polarizzazione 
                              double Phi_collins = Phi_HP_Rad + Phi_s; 
                              double Phi_sivers = Phi_HP_Rad - Phi_s;
                              if(Phi_collins< -TMath::Pi()){
                                Phi_collins+= 2*TMath::Pi();
                              }
                              else if(Phi_collins> TMath::Pi()){
                                Phi_collins-= 2*TMath::Pi();
                              }
                              double SinCollins = TMath::Sin(Phi_collins);
                              double F_polariz = PolarizFunction(Phi_collins, Phi_sivers);
                              double rho = RandomProb();
                              if(z_DA_pion <=  1){
                                  //xQplane->Fill(xbj_DA_pion, Q2_DA_pion);
                                  //pion_PhT->Fill(PhT_pion);
                                  //pion_PT->Fill(P_T);
                                  if(rho > F_polariz){
                                    spinDown += 1;
                                    TVector3 t1_1 = x_axis_SP.Cross(ProtonSpinDown);
                                    double t1 = t1_1 * z_axis_SP;
                                    double t2 = x_axis_SP.Dot(ProtonSpinDown);
                                    Phi_s_Down = std::atan2(t1, t2);
                                    Delta_Down = Phi_HP_Rad - Phi_s_Down;
                                    Delta_Down_Coll = Phi_HP_Rad + Phi_s_Down;
                                    if(Delta_Down < -TMath::Pi()){
                                      Delta_Down += 2*TMath::Pi();
                                    }
                                    else if(Delta_Down > TMath::Pi()){
                                      Delta_Down -= 2*TMath::Pi();
                                    }
                                    if(Delta_Down_Coll < -TMath::Pi()){
                                      Delta_Down_Coll += 2*TMath::Pi();
                                    }
                                    else if(Delta_Down_Coll > TMath::Pi()){
                                      Delta_Down_Coll -= 2*TMath::Pi();
                                    }
                                    SinPhiDown = TMath::Sin(Delta_Down);
                                    SinPhiDown_Coll = TMath::Sin(Delta_Down_Coll);
                                    DeltaPhi_Polarized->Fill(Delta_Down, SinPhiDown);
                                    DeltaPhi_Polarized_Area->Fill(Delta_Down);
                                    DeltaPhi_Pol_Collins->Fill(Delta_Down_Coll, SinPhiDown_Coll);
                                    DeltaPhi_Pol_CollinsArea->Fill(Delta_Down_Coll);
                                  }
                                  if(rho <= F_polariz){
                                    spinUp += 1;
                                    TVector3 t3_1 = x_axis_SP.Cross(ProtonSpin);
                                    double t3 = t3_1 * z_axis_SP;
                                    double t4 = x_axis_SP.Dot(ProtonSpin);
                                    Phi_s_Up = std::atan2(t3, t4);
                                    Delta_Up = Phi_HP_Rad - Phi_s_Up;
                                    Delta_Up_Coll = Phi_HP_Rad + Phi_s_Up;
                                    if(Delta_Up < -TMath::Pi()){
                                      Delta_Up += 2*TMath::Pi();
                                    }
                                    else if(Delta_Up > TMath::Pi()){
                                      Delta_Up -= 2*TMath::Pi();
                                    }
                                    if(Delta_Up_Coll < -TMath::Pi()){
                                      Delta_Up_Coll += 2*TMath::Pi();
                                    }
                                    else if(Delta_Up_Coll > TMath::Pi()){
                                      Delta_Up_Coll -= 2*TMath::Pi();
                                    }
                                    SinPhiUp = TMath::Sin(Delta_Up);
                                    SinPhiUp_Coll = TMath::Sin(Delta_Up_Coll);
                                    DeltaPhi_Polarized->Fill(Delta_Up, SinPhiUp);
                                    DeltaPhi_Polarized_Area->Fill(Delta_Up);
                                    DeltaPhi_Pol_Collins->Fill(Delta_Up_Coll, SinPhiUp_Coll);
                                    DeltaPhi_Pol_CollinsArea->Fill(Delta_Up_Coll);
                                  }
                                  pion_MomVsQ2->Fill(mom_pion, Q2_DA_pion);
                                  pion_MomVsEta->Fill(pionEta, mom_pion);
                                  pion_PhTvsQ2->Fill(PhT_pion, Q2_DA_pion);
                                  pion_PhTvsEta->Fill(pionEta, PhT_pion);
                                  pion_MomVsPhi->Fill(mom_pion, pionPhi);
                                  pion_PhTvsZ->Fill(z_DA_pion, PhT_pion);
                                  pion_MomVsTheta->Fill(pionTheta, mom_pion);
                                  pion_MomVsTheta_HP->Fill(mom_pion,Theta_HP_Deg);
                                  pion_MomVsPhi_HP->Fill(Phi_HP_Deg, mom_pion);
                                  pion_ThetaVsPhi->Fill(pionTheta, pionPhi);
                                  //pion_ThetaVsPhi2->Fill(pionTheta,pionPhi);
                                  pion_Theta_HPVsPhi_HP->Fill(Phi_HP_Deg, Theta_HP_Deg);
                                  //pion_Theta_HPVsPhi_HP2->Fill(Phi_HP_Deg, Theta_HP_Deg);
                                  pion_PhTvsPhi->Fill(PhT_pion, pionPhi);
                                  pion_PhTvsTheta->Fill(pionTheta, PhT_pion);
                                  pion_PhTvsPhi_HP->Fill(Phi_HP_Deg, PhT_pion);
                                  pion_PhTvsTheta_HP->Fill(PhT_pion,Theta_HP_Deg);
                                  //pion_MomVsPhiVsTheta->Fill(pionTheta, mom_pion, pionPhi);
                                  pion_PhTvsMom->Fill(PhT_pion,mom_pion);
                                  pion_PhTvsMom2->Fill(PhT_pion,mom_pion);
                                  pion_ZvsQ2->Fill(z_DA_pion, Q2_DA_pion);
                                  pion_ZvsMom->Fill(z_DA_pion, mom_pion);
                                  pion_ZvsPhi->Fill(z_DA_pion, pionPhi);
                                  pion_ZvsTheta->Fill(pionTheta, z_DA_pion);
                                  pion_ZvsPhi_HP->Fill(Phi_HP_Deg,z_DA_pion);
                                  pion_ZvsTheta_HP->Fill(z_DA_pion,Theta_HP_Deg);
                                  pion_XbVsMom->Fill(xbj_DA_pion, mom_pion);
                                  pion_XbVsPhT->Fill(xbj_DA_pion, PhT_pion);
                                  pion_XbVsZ->Fill(xbj_DA_pion, z_DA_pion);
                                  pion_XbVsPhi->Fill(xbj_DA_pion, pionPhi);
                                  //pion_XbVsPhi2->Fill(pionPhi, xbj_DA_pion);
                                  pion_XbVsTheta->Fill(pionTheta, xbj_DA_pion);
                                  pion_XbVsPhi_HP->Fill(Phi_HP_Deg, xbj_DA_pion);
                                  pion_XbVsTheta_HP->Fill(xbj_DA_pion, Theta_HP_Deg);
                                  //pion_XbVsTheta_HP2->Fill(Theta_HP_Deg, xbj_DA_pion);
                                  pion_XbVsQ2->Fill(xbj_DA_pion, Q2_DA_pion);
                                  //pion_PhTvsZvsXb->Fill(xbj_DA_pion, z_DA_pion, PhT_pion);
                                  pion_XbVsEta->Fill(pionEta, xbj_DA_pion);
                                  //z_axisSP->Fill(z_axis_z, z_axis_x, z_axis_y);
                                  pion_P_TvsPhi->Fill(P_T, pionPhi);
                                  pion_P_TvsTheta->Fill(pionTheta, P_T);
                                  pion_P_TvsPhi_HP->Fill(Phi_HP_Deg, P_T);
                                  pion_P_TvsTheta_HP->Fill(P_T, Theta_HP_Deg);
                                  pion_P_TvsEta->Fill(pionEta, P_T);
                                  pion_P_TvsMom->Fill(P_T,mom_pion);
                                  pion_P_TvsMom2->Fill(P_T,mom_pion);
                                  pion_P_TvsQ2->Fill(P_T, Q2_DA_pion);
                                  pion_P_TvsZ->Fill(z_DA_pion, P_T);
                                  pion_P_TvsXb->Fill(xbj_DA_pion, P_T);
                                  pion_ThetaVsTheta_HP->Fill(Theta_HP_Deg, pionPhi);
                                  pion_PhiVsPhi_HP->Fill(pionTheta, Phi_HP_Deg);
                                  pion_Q2vsEta->Fill(pionEta, Q2_DA_pion);
                                  //axisHP->Fill(z, x, y);
                                  pion_EtaVsPhi_HP->Fill(Phi_HP_Deg, pionEta);
                                  pion_ThetaVsPhi_HP->Fill(Phi_HP_Deg, pionPhi);
                                  Phi_ss->Fill(Phi_s_Deg);
                                  Phi_ss2->Fill(phi_s22_Deg);
                                  Phi_h->Fill(Phi_HP_Deg);
                                  Phi_h2->Fill(phi_hh_Deg);
                                  pion_Phi_HPvsPhi_s->Fill(Phi_HP_Deg, Phi_s_Deg);
                                  pion_Phi_HPvsPhi_s2->Fill(Phi_acca_Deg, Phi_s_Deg);
                                  SinPhi_plot->Fill(SinPhi);
                                  Cos2Phi_plot->Fill(Cos2Phi_h);
                                  SinPhivsDelta->Fill(DeltaPhi, SinPhi);
                                  Cos2PhivsPhi->Fill(2*Phi_HP_Rad, Cos2Phi_h);
                                  pion_SinPhiVsZ->Fill(SinPhi, z_DA_pion);
                                  pion_SinPhiVsPhT->Fill(SinPhi, PhT_pion);
                                  pion_SinPhiVsMom->Fill(SinPhi, mom_pion);
                                  pion_SinPhiVsXb->Fill(SinPhi, xbj_DA_pion);
                                  pion_SinPhiVsQ2->Fill(SinPhi,Q2_DA_pion);
                                  pion_DeltaPhiVsMom->Fill(DeltaPhi_Rad, mom_pion);
                                  pion_DeltaPhiVsPhT->Fill(DeltaPhi_Rad, PhT_pion);
                                  pion_DeltaPhiVsTheta->Fill(DeltaPhi_Rad,pionPhi);
                                  pion_DeltaPhiVsXb->Fill(DeltaPhi_Rad, xbj_DA_pion);
                                  pion_DeltaPhiVsZ->Fill(DeltaPhi_Rad, z_DA_pion);
                                  z_plot->Fill(z_DA_pion);
                                  DeltaPhi_plot->Fill(DeltaPhi);
                                  DeltaPhi_plot2->Fill(DeltaPhi2);
                                  DeltaPhi_plot3->Fill(DeltaPhi, SinPhi);
                                  DeltaPhi_sinArea->Fill(DeltaPhi);
                                  P_over_Q->Fill(Rapporto);
                                  DeltaPhi_Collins->Fill(Phi_collins, SinCollins);
                                  DeltaPhi_CollinsArea->Fill(Phi_collins);
                                  //RealdRICH_xy->Fill(dRICH_x,dRICH_y);
                                  if(xbj_DA_pion <= 0.1 && xbj_DA_pion >= 0.01){
                                    DeltaPhi_xLow->Fill(DeltaPhi, SinPhi);
                                  }
                                  if(xbj_DA_pion >= 0.1){
                                    DeltaPhi_xHigh->Fill(DeltaPhi, SinPhi);
                                  }
                                  if(pion_EpZ >= 1980){
                                    dRICH_xy->Fill(x_drich, y_drich);
                                  }
                                  if(pdg==211){
                                    pion_Q2_dist->Fill(Q2_DA_pion);
                                    pion_X_dist->Fill(xbj_DA_pion);
                                    pion_z_dist->Fill(z_DA_pion);
                                    pion_PhT_dist->Fill(PhT_pion);
                                  }
                              }
                            }
                        }
                    }
                  }
                }  
              
              
            }
        }
    } 
    

    // CALCOLO COEFFICIENTI DI FOURIER
    double A, a0, a1, b1;
    int fBins = DeltaPhi_plot->GetNbinsX();
    double f_deltaPhi = DeltaPhi_plot->GetBinWidth(1);
    double normFactor = f_deltaPhi / TMath::Pi();
    for (int i = 1; i <= fBins; ++i) {
        double phi_f = DeltaPhi_plot->GetBinCenter(i); // Centro del bin
        double f_phi = DeltaPhi_plot->GetBinContent(i); // Conteggio nel bin
        a0 += (f_phi * normFactor);
        a1 += f_phi * TMath::Cos(phi_f) * normFactor;
        b1 += f_phi * TMath::Sin(phi_f) * normFactor;
    }
    TCanvas* fourier = new TCanvas("Fourier", "Analisi di Fourier", 800, 600);
    DeltaPhi_plot->Draw("E");
    /*
    TF1* fitFourier = new TF1("fitFourier", "0.25*[0] + [1]*cos(x) + [2]*sin(x)", -2*TMath::Pi(), 2*TMath::Pi());
    fitFourier->SetParameters(a0, a1, b1); 
    fitFourier->SetLineColor(kRed);
    fitFourier->Draw("SAME");
    */
    fourier->Update();
    fourier->Write();
    delete fourier;
    TCanvas* fourier2 = new TCanvas("Fourier2", "Analisi di Fourier", 800, 600);
    DeltaPhi_plot2->Draw("E");
    fourier2->Update();
    fourier2->Write();
    delete fourier2;

    // canvas per siver
    TCanvas* sivers = new TCanvas("Sivers", "Sivers effect analysis", 800, 600);
    DeltaPhi_plot3->Draw("E");
    sivers->Update();
    sivers->Write();
    delete sivers;

    TCanvas* siversPol = new TCanvas("Sivers_Polarized", "Sivers effect analysis", 800, 600);
    DeltaPhi_Polarized->Draw("E");
    siversPol->Update();
    siversPol->Write();
    delete siversPol;


    /*
    // PER IL GRAFICO IN 3D
    //gStyle->SetPalette(kRainBow);  // Usa la palette dei colori "RainBow"
    gStyle->SetOptStat(0);         // Disattiva le informazioni statistiche
    // Creazione di una nuova canvas
    TCanvas *c1 = new TCanvas("Particle_XvsZvsPhT", "3D Histogram", 800, 600);
    c1->SetLogx();  
    c1->SetLogy();  
    c1->SetLogz();  
    // Disegna l'istogramma 3D con i colori sui bin
    pion_PhTvsZvsXb->Draw("BOX2Z");
    // Aggiorna il canvas per visualizzare tutto correttamente
    c1->Update();
    c1->Write();
    TCanvas *c2 = new TCanvas("Particle_ThetavsMomvsPhi", "3D Histogram", 800, 600);
    //c1->SetLogx();  
    c2->SetLogy();  
    //c1->SetLogz();  
    // Disegna l'istogramma 3D con i colori sui bin
    pion_MomVsPhiVsTheta->Draw("BOX2Z");
    // Aggiorna il canvas per visualizzare tutto correttamente
    c2->Update();
    c2->Write();
    */

    /*
    TCanvas *c3 = new TCanvas("beam axis", "Grafico in coordinate cartesiane", 800, 600);
    graphProton->SetMarkerColor(kBlue);
    graphElectron->SetMarkerColor(kRed);
    graphProton->SetTitle("Protoni - Coordinate cartesiane; X; Y; Z");
    graphElectron->SetTitle("Elettroni - Coordinate cartesiane; X; Y; Z");

    c3->Divide(2, 1);
    c3->cd(1);
    graphProton->Draw("P0 COLZ");  // Grafico per i protoni

    c3->cd(2);
    graphElectron->Draw("P0 COLZ");  // Grafico per gli elettroni
    c3->Update();
    c3->Write();
    */
    TCanvas *Mom = new TCanvas("MomVsAngles", "Correlation of the hadron momentum", 1200, 800);
    Mom->Divide(2,2);
    Mom->cd(1);
    pion_MomVsPhi->Draw();
    TPaveStats *stats1 = (TPaveStats*)pion_MomVsPhi->GetListOfFunctions()->FindObject("stats");
    if (stats1) {
        stats1->SetX1NDC(0.65); // Posizione X
        stats1->SetX2NDC(0.8); // Posizione X
        stats1->SetY1NDC(0.65); // Posizione Y
        stats1->SetY2NDC(0.8); // Posizione Y
        stats1->Draw(); // Ridisegna il pannello delle statistiche
    }
    Mom->cd(4);
    pion_MomVsPhi_HP->Draw();
    Mom->cd(2);
    pion_MomVsTheta->Draw();
    Mom->cd(3);
    pion_MomVsTheta_HP->Draw();

    Mom->Update();
    Mom->Write();
    delete Mom;
    
    // CANVAS DEI GRAFICI LOGARITMICI
    gStyle->SetOptStat(1111);
    TCanvas *c3 = new TCanvas("Q2vsMom", "Mom vs Q^2", 800, 600);
    c3->SetLogy();
    pion_MomVsQ2->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats = (TPaveStats*)pion_MomVsQ2->FindObject("stats");
    if (stats) {  
        stats->SetX1NDC(0.72);  // Sposta il lato sinistro al 70% della larghezza
        stats->SetX2NDC(0.88);  // Sposta il lato destro al 90% della larghezza
        stats->SetY1NDC(0.72);  // Sposta il lato inferiore al 20% dell'altezza
        stats->SetY2NDC(0.88);  // Sposta il lato superiore al 40% dell'altezza
    }
    gPad->Modified();
    gPad->Update();
    c3->Update();
    c3->Write();
    delete c3;

    TCanvas *c4 = new TCanvas("Q2vsPhT", "Z vs Q^2", 800, 600);
    c4->SetLogy();
    pion_PhTvsQ2->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats2 = (TPaveStats*)pion_PhTvsQ2->FindObject("stats");
    if (stats2) {  
        stats2->SetX1NDC(0.72);  // Sposta il lato sinistro al 70% della larghezza
        stats2->SetX2NDC(0.88);  // Sposta il lato destro al 90% della larghezza
        stats2->SetY1NDC(0.72);  // Sposta il lato inferiore al 20% dell'altezza
        stats2->SetY2NDC(0.88);  // Sposta il lato superiore al 40% dell'altezza
    }
    gPad->Modified();
    gPad->Update();
    c4->Update();
    c4->Write();
    delete c4;

    TCanvas *c14 = new TCanvas("Q2vsP_T", "Z vs Q^2", 800, 600);
    c14->SetLogy();
    pion_P_TvsQ2->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats3 = (TPaveStats*)pion_P_TvsQ2->FindObject("stats");
    if (stats3) {  
        stats3->SetX1NDC(0.72);  // Sposta il lato sinistro al 70% della larghezza
        stats3->SetX2NDC(0.88);  // Sposta il lato destro al 90% della larghezza
        stats3->SetY1NDC(0.72);  // Sposta il lato inferiore al 20% dell'altezza
        stats3->SetY2NDC(0.88);  // Sposta il lato superiore al 40% dell'altezza
    }
    gPad->Modified();
    gPad->Update();
    c14->Update();
    c14->Write();
    delete c14;

    TCanvas *c5 = new TCanvas("Q2vsZ", "Z vs Q^2", 800, 600);
    c5->SetLogy();
    pion_ZvsQ2->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats4 = (TPaveStats*)pion_ZvsQ2->FindObject("stats");
    if (stats4) {  
        stats4->SetX1NDC(0.72);  // Sposta il lato sinistro al 70% della larghezza
        stats4->SetX2NDC(0.88);  // Sposta il lato destro al 90% della larghezza
        stats4->SetY1NDC(0.72);  // Sposta il lato inferiore al 20% dell'altezza
        stats4->SetY2NDC(0.88);  // Sposta il lato superiore al 40% dell'altezza
    }
    gPad->Modified();
    gPad->Update();
    c5->Update();
    c5->Write();
    delete c5;

    TCanvas *c6 = new TCanvas("Q2vsXb", "Z vs Q^2", 800, 600);
    c6->SetLogy();
    c6->SetLogx();
    pion_XbVsQ2->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats5 = (TPaveStats*)pion_XbVsQ2->FindObject("stats");
    if (stats5) {  
        stats5->SetX1NDC(0.12);  // Sposta il lato sinistro al 70% della larghezza
        stats5->SetX2NDC(0.28);  // Sposta il lato destro al 90% della larghezza
        stats5->SetY1NDC(0.72);  // Sposta il lato inferiore al 20% dell'altezza
        stats5->SetY2NDC(0.88);  // Sposta il lato superiore al 40% dell'altezza
    }
    gPad->Modified();
    gPad->Update();
    c6->Update();
    c6->Write();
    delete c6;

    TCanvas *c17 = new TCanvas("Q2vsEta", "Z vs Q^2", 800, 600);
    c17->SetLogy();
    pion_Q2vsEta->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats6 = (TPaveStats*)pion_Q2vsEta->FindObject("stats");
    if (stats6) {  
        stats6->SetX1NDC(0.12);  // Sposta il lato sinistro al 70% della larghezza
        stats6->SetX2NDC(0.28);  // Sposta il lato destro al 90% della larghezza
        stats6->SetY1NDC(0.72);  // Sposta il lato inferiore al 20% dell'altezza
        stats6->SetY2NDC(0.88);  // Sposta il lato superiore al 40% dell'altezza
    }
    gPad->Modified();
    gPad->Update();
    c17->Update();
    c17->Write();
    delete c17;

    TCanvas *c7 = new TCanvas("XbVsMom", "Z vs Q^2", 800, 600);
    c7->SetLogx();
    pion_XbVsMom->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats7 = (TPaveStats*)pion_XbVsMom->FindObject("stats");
    if (stats7) {  
        stats7->SetX1NDC(0.12);  // Sposta il lato sinistro al 70% della larghezza
        stats7->SetX2NDC(0.28);  // Sposta il lato destro al 90% della larghezza
        stats7->SetY1NDC(0.72);  // Sposta il lato inferiore al 20% dell'altezza
        stats7->SetY2NDC(0.88);  // Sposta il lato superiore al 40% dell'altezza
    }
    gPad->Modified();
    gPad->Update();
    c7->Update();
    c7->Write();
    delete c7;

    TCanvas *c8 = new TCanvas("XbVsPhT", "Z vs Q^2", 800, 600);
    c8->SetLogx();
    pion_XbVsPhT->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats8 = (TPaveStats*)pion_XbVsPhT->FindObject("stats");
    if (stats8) {  
        stats8->SetX1NDC(0.12);  // Sposta il lato sinistro al 70% della larghezza
        stats8->SetX2NDC(0.28);  // Sposta il lato destro al 90% della larghezza
        stats8->SetY1NDC(0.72);  // Sposta il lato inferiore al 20% dell'altezza
        stats8->SetY2NDC(0.88);  // Sposta il lato superiore al 40% dell'altezza
    }
    gPad->Modified();
    gPad->Update();
    c8->Update();
    c8->Write();
    delete c8;

    TCanvas *c15 = new TCanvas("XbVsP_T", "Z vs Q^2", 800, 600);
    c15->SetLogx();
    pion_P_TvsXb->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats15 = (TPaveStats*)pion_P_TvsXb->FindObject("stats");
    if (stats15) {  
        stats15->SetX1NDC(0.72);  // Sposta il lato sinistro al 70% della larghezza
        stats15->SetX2NDC(0.88);  // Sposta il lato destro al 90% della larghezza
        stats15->SetY1NDC(0.72);  // Sposta il lato inferiore al 20% dell'altezza
        stats15->SetY2NDC(0.88);  // Sposta il lato superiore al 40% dell'altezza
    }
    gPad->Modified();
    gPad->Update();
    c15->Update();
    c15->Write();
    delete c15;

    TCanvas *c9 = new TCanvas("XbVsZ", "Z vs Q^2", 800, 600);
    c9->SetLogx();
    pion_XbVsZ->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats9 = (TPaveStats*)pion_XbVsZ->FindObject("stats");
    if (stats9) {  
        stats9->SetX1NDC(0.72);  // Sposta il lato sinistro al 70% della larghezza
        stats9->SetX2NDC(0.88);  // Sposta il lato destro al 90% della larghezza
        stats9->SetY1NDC(0.72);  // Sposta il lato inferiore al 20% dell'altezza
        stats9->SetY2NDC(0.88);  // Sposta il lato superiore al 40% dell'altezza
    }
    gPad->Modified();
    gPad->Update();
    c9->Update();
    c9->Write();
    delete c9;
    
    TCanvas *c16 = new TCanvas("XbVsEta", "Z vs Q^2", 800, 600);
    c16->SetLogy();
    pion_XbVsEta->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats16 = (TPaveStats*)pion_XbVsEta->FindObject("stats");
    if (stats16) {  
        stats16->SetX1NDC(0.12);  // Sposta il lato sinistro al 70% della larghezza
        stats16->SetX2NDC(0.28);  // Sposta il lato destro al 90% della larghezza
        stats16->SetY1NDC(0.72);  // Sposta il lato inferiore al 20% dell'altezza
        stats16->SetY2NDC(0.88);  // Sposta il lato superiore al 40% dell'altezza
    }
    gPad->Modified();
    gPad->Update();
    c16->Update();
    c16->Write();
    delete c16;

    TCanvas *c10 = new TCanvas("XbVsTheta", "Z vs Q^2", 800, 600);
    c10->SetLogx();
    pion_XbVsPhi->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats10 = (TPaveStats*)pion_XbVsPhi->FindObject("stats");
    if (stats10) {  
        stats10->SetX1NDC(0.72);  // Sposta il lato sinistro al 70% della larghezza
        stats10->SetX2NDC(0.88);  // Sposta il lato destro al 90% della larghezza
        stats10->SetY1NDC(0.72);  // Sposta il lato inferiore al 20% dell'altezza
        stats10->SetY2NDC(0.88);  // Sposta il lato superiore al 40% dell'altezza
    }
    gPad->Modified();
    gPad->Update();
    c10->Update();
    c10->Write();
    delete c10;

    TCanvas *c11 = new TCanvas("XbVsPhi", "Z vs Q^2", 800, 600);
    c11->SetLogy();
    pion_XbVsTheta->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats11 = (TPaveStats*)pion_XbVsTheta->FindObject("stats");
    if (stats11) {  
        stats11->SetX1NDC(0.72);  // Sposta il lato sinistro al 70% della larghezza
        stats11->SetX2NDC(0.88);  // Sposta il lato destro al 90% della larghezza
        stats11->SetY1NDC(0.72);  // Sposta il lato inferiore al 20% dell'altezza
        stats11->SetY2NDC(0.88);  // Sposta il lato superiore al 40% dell'altezza
    }
    gPad->Modified();
    gPad->Update();
    c11->Update();
    c11->Write();
    delete c11;

    TCanvas *c12 = new TCanvas("XbVsTheta_HP", "Z vs Q^2", 800, 600);
    c12->SetLogx();
    pion_XbVsTheta_HP->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats12 = (TPaveStats*)pion_XbVsTheta_HP->FindObject("stats");
    if (stats12) {  
        stats12->SetX1NDC(0.72);  // Sposta il lato sinistro al 70% della larghezza
        stats12->SetX2NDC(0.88);  // Sposta il lato destro al 90% della larghezza
        stats12->SetY1NDC(0.72);  // Sposta il lato inferiore al 20% dell'altezza
        stats12->SetY2NDC(0.88);  // Sposta il lato superiore al 40% dell'altezza
    }
    gPad->Modified();
    gPad->Update();
    c12->Update();
    c12->Write();
    delete c12;

    TCanvas *c13 = new TCanvas("XbVsPhi_HP", "Z vs Q^2", 800, 600);
    c13->SetLogy();
    pion_XbVsPhi_HP->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats13 = (TPaveStats*)pion_XbVsPhi_HP->FindObject("stats");
    if (stats13) {  
        stats13->SetX1NDC(0.72);  
        stats13->SetX2NDC(0.88);  
        stats13->SetY1NDC(0.72);  
        stats13->SetY2NDC(0.88);  
    }
    gPad->Modified();
    gPad->Update();
    c13->Update();
    c13->Write();
    delete c13;
    /*
    if (!RealdRICH_xy) {
    std::cerr << "Error: RealdRICH_xy is not initialized!" << std::endl;
    return;
    }
    */
    /*
    TCanvas *c18 = new TCanvas("dRICH_xy", "dRICH xy Plane", 800, 600);
    dRICH_xy->Draw("COLZ");
    gPad->Update(); 
    TPaveStats *stats18 = (TPaveStats*)dRICH_xy->FindObject("stats");
    if (stats18) {  
        stats18->SetX1NDC(0.72);  
        stats18->SetX2NDC(0.88);  
        stats18->SetY1NDC(0.72);  
        stats18->SetY2NDC(0.88);
    }
    TLine *lineX = new TLine(-300, 0, 300, 0);  // Linea orizzontale (y = 0)
    TLine *lineY = new TLine(0, -300, 0, 300);  // Linea verticale (x = 0)
    // Imposta lo stile tratteggiato
    lineX->SetLineStyle(2);  // Stile tratteggiato
    lineY->SetLineStyle(2);  // Stile tratteggiato
    // Disegna le linee sopra il TH2D
    lineX->Draw("SAME");
    lineY->Draw("SAME");
    gPad->Modified();
    //gPad->Update();
    c18->Update();
    c18->Write();

    TCanvas *c19 = new TCanvas("dRICH_xy_Rec", "dRICH xy Plane (Reco)", 800, 600);
    dRICH_xy_Rec->Draw("COLZ");
    TLine *lineX2 = new TLine(-300, 0, 300, 0);  // Linea orizzontale (y = 0)
    TLine *lineY2 = new TLine(0, -300, 0, 300);  // Linea verticale (x = 0)
    // Imposta lo stile tratteggiato
    lineX2->SetLineStyle(2);  // Stile tratteggiato
    lineY2->SetLineStyle(2);  // Stile tratteggiato
    // Disegna le linee sopra il TH2D
    lineX2->Draw("SAME");
    lineY2->Draw("SAME");
    gPad->Update(); 
    TPaveStats *stats19 = (TPaveStats*)dRICH_xy_Rec->FindObject("stats");
    if (stats19) {  
        stats19->SetX1NDC(0.72);  
        stats19->SetX2NDC(0.88);  
        stats19->SetY1NDC(0.72);  
        stats19->SetY2NDC(0.88);
    }
    gPad->Modified();
    //gPad->Update();
    c19->Update();
    c19->Write();

    TCanvas *c20 = new TCanvas("dRICH_xy_PID", "dRICH xy Plane (Reco x PID)", 800, 600);
    dRICH_xy_PID->Draw("COLZ");
    TLine *lineX3 = new TLine(-300, 0, 300, 0);  // Linea orizzontale (y = 0)
    TLine *lineY3 = new TLine(0, -300, 0, 300);  // Linea verticale (x = 0)
    // Imposta lo stile tratteggiato
    lineX3->SetLineStyle(2);  // Stile tratteggiato
    lineY3->SetLineStyle(2);  // Stile tratteggiato
    // Disegna le linee sopra il TH2D
    lineX3->Draw("SAME");
    lineY3->Draw("SAME");
    gPad->Update(); 
    TPaveStats *stats20 = (TPaveStats*)dRICH_xy_PID->FindObject("stats");
    if (stats20) {  
        stats20->SetX1NDC(0.72);  
        stats20->SetX2NDC(0.88);  
        stats20->SetY1NDC(0.72);  
        stats20->SetY2NDC(0.88);
    }
    gPad->Modified();
    //gPad->Update();
    c20->Update();
    c20->Write();
    */
    
    //std::cout << "particelle generate: " << count << std::endl;
    //std::cout << "particelle generate nel range di rapidita': " << countEta << std::endl;
    std::cout << "the acceptance is: " << countEta / count  << std::endl;
    double A_ut = (spinUp-spinDown)/(spinUp+spinDown);
    std::cout << "SpinUp: " << spinUp << std::endl;
    std::cout << "SpinDown: " << spinDown << std::endl;
    std::cout << "A_UT = " << A_ut << std::endl;
    std::cout << "______________________________________________________________________________________" << std::endl;
    std::cout << " " << std::endl;
    //std::cout << "dRICHx: " << dRICHx << ", dRICHy: " << dRICHy << std::endl;
    

    /*
    std::cout << "Coefficiente A (assimmetria): " << b1/(0.25*a0) << std::endl;
    std::cout << "Coefficiente a0 (termine medio): " << 0.25*a0 << std::endl;
    std::cout << "Coefficiente a1 (cos): " << a1 << std::endl;
    std::cout << "Coefficiente b1 (sin): " << b1 << std::endl;
    */
    //std::cout << "ProtonBeam = (" << ProtonBeam.Px() << "," << ProtonBeam.Py() << "," << ProtonBeam.Pz() << "," << ProtonBeam.E() << ")" << std::endl;
    //std::cout << "elettroni lanciati: " << count_el << std::endl;
    /*
    for (const auto& status : uniqueStatuses) {
        std::cout << "lo status e': " << status << std::endl;
    }
    
    for (const auto& index : uniqueParentsIndex) {
        std::cout << "l'indice del genitore e': " << index << std::endl;
    }
    */

    // _____________________________________________________________________________________________________

    ofile->Write(); // Write histograms to file
    ofile->Close(); // Close output file

    mychain->Delete();
    ofile->Delete();
  }


int main11() {

  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1808.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1840.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1909.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out11.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main12() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1930.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1929.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1928.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out12.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main13() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1926.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1925.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1928.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out13.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main14() {
  // il 1924 non esiste
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1923.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1923.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1922.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out14.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main15() {
  // anche il 20 manca
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1921.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1921.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1919.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out15.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main16() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1918.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1917.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1916.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out16.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main17() {
  // manca il 13
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1915.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1914.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1914.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out17.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main18() {
  // 11
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1912.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1910.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1910.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out18.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main19() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1909.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1908.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1907.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out19.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main20() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1906.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1905.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1904.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out20.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main21() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1903.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1902.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1901.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out21.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main22() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1900.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1899.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1898.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out22.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main23() {
  // 96
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1898.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1897.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1897.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out23.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main24() {
  // 93
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1895.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1894.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1894.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out24.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main25() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1892.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1891.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1890.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out25.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main26() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1889.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1888.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1886.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out26.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main27() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1886.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1885.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1884.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out27.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main28() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1883.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1882.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1881.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out28.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main29() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1880.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1879.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1878.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out29.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main30() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1877.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1876.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1875.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out30.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main31() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1874.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1873.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1872.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out31.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main32() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1871.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1870.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1869.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out32.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main33() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1868.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1867.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1866.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out33.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main34() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1865.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1864.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1863.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out34.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main35() {
  // 60
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1862.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1861.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1861.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out35.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main36() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1859.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1858.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1857.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out36.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main37() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1856.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1855.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1854.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out37.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main38() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1853.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1852.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1851.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out38.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main39() {
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1850.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1849.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1848.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out39.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

int main40() {
  // 47 e 46
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1845.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1844.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1843.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out40.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);

  return 0;

}

// da 37 a 32 non esistono, 38 e 31 si
int main41() {
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1842.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1841.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1840.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out41.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);
  return 0;
}

int main42() {
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1839.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1838.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1831.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out42.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);
  return 0;
}

int main43() {
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1830.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1829.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1828.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out43.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);
  return 0;
}

int main44() {
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1824.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1826.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1825.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out44.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);
  return 0;
}

int main45() {
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1824.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1823.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1822.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out45.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);
  return 0;
}

int main46() {
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1810.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1820.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1819.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out46.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);
  return 0;
}
// no: 11, 13, 21, 27
int main47() {
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1818.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1817.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1816.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out47.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);
  return 0;
}
int main48() {
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1815.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1814.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1812.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out48.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);
  return 0;
}
// 94, 95, 
int main49() {
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1810.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1809.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1808.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out49.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);
  return 0;
}

int main50() {
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1807.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1806.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1805.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out50.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);
  return 0;
}

int main51() {
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1803.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1804.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1802.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out51.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);
  return 0;
}

int main52() {
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1801.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1800.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1799.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out52.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);
  return 0;
}

int main53() {
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1798.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1797.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1796.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out53.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);
  return 0;
}

int main54() {
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1793.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1792.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1798.eicrecon.tree.edm4eic.root";
  const char* outputFile = "out54.histAbs.root";
  pino(inputFile1, inputFile2, inputFile3, outputFile);
  return 0;
}
