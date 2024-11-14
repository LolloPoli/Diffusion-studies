#include <TH1D.h>
#include <TH2F.h>
#include <TGraphErrors.h>
//#include <LHAPDF/LHAPDF.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TAxis.h>
#include <vector>
#include <TMath.h>
#include "TLorentzVector.h"
#include <cmath>
// per far funzionare il programma
// . /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-centos7-gcc11-opt/setup.sh
// root
// .L taskB.C
// .L run_all.C
// run()
// NON SO MAGARI SERVONO PER LA NNPDF
// gSystem->AddIncludePath("/cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-centos7-gcc11-opt/include");
//gSystem->AddLinkPath("/cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-centos7-gcc11-opt/lib");

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
/*
// 1889, 1888, 1887, 1886, 1885, 1884, 1883, 1882, 1881, 1880, 1822, 1803, 1802, 1801, 1800
// pythia8CCDIS_18x275_minQ2=1000_beamEffects_xAngle=-0.025_hiDiv_5.2026.eicrecon.tree.edm4eic.root
// pythia_ep_noradcor_10x275_q2_0.000000001_1.0_run15.ab.0199.eicrecon.tree.edm4eic.root
// pythia_ep_noradcor_18x275_q2_0.000000001_1.0_run9.ab.1887.eicrecon.tree.edm4eic.root
void pino(TString infile1="pythia_ep_noradcor_18x275_q2_0.000000001_1.0_run9.ab.1887.eicrecon.tree.edm4eic.root", 
TString infile2="pythia_ep_noradcor_18x275_q2_0.000000001_1.0_run9.ab.1888.eicrecon.tree.edm4eic.root",
TString infile3="pythia_ep_noradcor_18x275_q2_0.000000001_1.0_run9.ab.1889.eicrecon.tree.edm4eic.root")
{ 
*/
/*
void pino(TString infile1="pythia_ep_noradcor_10x275_q2_0.000000001_1.0_run15.ab.0199.eicrecon.tree.edm4eic.root", 
TString infile2="pythia_ep_noradcor_10x275_q2_0.000000001_1.0_run15.ab.0198.eicrecon.tree.edm4eic.root",
TString infile3="pythia_ep_noradcor_10x275_q2_0.000000001_1.0_run15.ab.0200.eicrecon.tree.edm4eic.root")
{
*/
// pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_4.1467.eicrecon.tree.edm4eic.root
// pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1688.eicrecon.tree.edm4eic.root
// da 1449 a 1469
void pino(const char* inputFile1, const char* inputFile2, const char* inputFile3, const char* outputFile){
/*
void pino(TString infile1="pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_4.1469.eicrecon.tree.edm4eic.root",
TString infile2="pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_4.1468.eicrecon.tree.edm4eic.root",
TString infile3="pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_4.1467.eicrecon.tree.edm4eic.root"){
*/

    // Set output file for the histograms 11-16
    TFile *ofile = TFile::Open(outputFile, "RECREATE");

    // Set up input file chain
    TChain *mychain = new TChain("events");
    mychain->Add(inputFile1);
    mychain->Add(inputFile2);
    mychain->Add(inputFile3);

    // Initialize reader
    TTreeReader tree_reader(mychain);

    // Get Particle Information
    TTreeReaderArray<int> partGenStat(tree_reader, "MCParticles.generatorStatus");
    TTreeReaderArray<float> partMomX(tree_reader, "MCParticles.momentum.x");
    TTreeReaderArray<float> partMomY(tree_reader, "MCParticles.momentum.y");
    TTreeReaderArray<float> partMomZ(tree_reader, "MCParticles.momentum.z");
    TTreeReaderArray<int> partPdg(tree_reader, "MCParticles.PDG");
    TTreeReaderArray<int> parentsIndex(tree_reader, "_MCParticles_parents.index");
    TTreeReaderArray<int> daughterIndex(tree_reader, "_MCParticles_daughters.index");
    TTreeReaderArray<unsigned int> par(tree_reader, "MCParticles.parents_end");


    // Get Reconstructed Track Information
    TTreeReaderArray<float> trackMomX(tree_reader, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float> trackMomY(tree_reader, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float> trackMomZ(tree_reader, "ReconstructedChargedParticles.momentum.z");
    TTreeReaderArray<int> recPdg(tree_reader, "ReconstructedChargedParticles.PDG");

    // Get Associations Between MCParticles and ReconstructedChargedParticles
    TTreeReaderArray<unsigned int> recoAssoc(tree_reader, "ReconstructedChargedParticleAssociations.recID");
    TTreeReaderArray<unsigned int> simuAssoc(tree_reader, "ReconstructedChargedParticleAssociations.simID");


    // SCATTERED ELECTRON
    TH1D *scatEl = new TH1D("scatEl", "scattered electron Mom; GeV^2", 80, 0, 20.);
    TH1D *recScatEl = new TH1D("RecScatElectron", "reconstruction of scattered electron Mom; GeV^2",80, 0, 18.);
    TH1D *scatElAngle = new TH1D("ScatElAngle", "Angle of the scattered electron; Theta", 80, 0., 180.);
    // GRAFICI DI ETA
    TH1D *truepionEta = new TH1D("truepionEta","Eta of charged Pions;Eta",80,-3.5,3.5);
    TH1D *RecpionEta = new TH1D("RecpionEta","Eta of reconstructed Pions+;Eta",80,-3.5,3.5);
    TH1D *truekaonEta = new TH1D("truekaonEta","Eta of charged Kaons;Eta",80,-3.5,3.5);
    TH1D *ReckaonEta = new TH1D("ReckaonEta","Eta of reconstructed K+;Eta",80,-3.5,3.5);
    TH1D *trueprotonEta = new TH1D("trueprotonEta","Eta of charged Protons;Eta",80,-3.5,3.5);
    TH1D *RecprotonEta = new TH1D("RecprotonEta","Eta of reconstructed Protons;Eta",80,-3.5,3.5);
    // GRAFICI DI PHI
    TH1D *truepionPhi = new TH1D("truepionPhi","Phi of charged Pions; Phi",80,0,25);
    TH1D *RecpionPhi = new TH1D("RecpionPhi","Phi of reconstructed Pions+; Phi",80,0,25);    
    TH1D *truekaonPhi = new TH1D("truekaonPhi","Phi of charged Kaons; Phi",80,0,25);
    TH1D *ReckaonPhi = new TH1D("ReckaonPhi","Phi of reconstructed Kaons+; Phi",80,0,25);  
    TH1D *trueprotonPhi = new TH1D("trueprotonPhi","Phi of charged Protons; Phi",80,0,25);
    TH1D *RecprotonPhi = new TH1D("RecprotonPhi","Phi of reconstructed Protons; Phi",80,0,25); 
    TH2D *PIDgen2rec = new TH2D("controlloPDG", "Reconstruction table; PDG gen; PDG rec",5, 0., 5., 5, 0., 5.);
    char part[6][3] = {"","e", "pi", "K", "p", "x"};
    for(int i = 1; i<6; i++){
      PIDgen2rec->GetXaxis()->SetBinLabel(i, part[i]);
      PIDgen2rec->GetYaxis()->SetBinLabel(i, part[i]);
    }
    // GRAFICI DI Q^2
    int nbins = 80;
    int nbon = 60;
    int nben = 60;
    int nbun = 40;
    double xmin_xbj = 5e-3;
    double xmax_xbj = 1;
    double xmin_Q2 = 0.9;
    double xmax_Q2 = 100.;
    double qm = 1;
    double qM = 1e4;
    double xM = 1;
    std::vector<double> log_bins_Q2 = CreateLogBinning(nbins, xmin_Q2, xmax_Q2);
    std::vector<double> log_bins_xbj = CreateLogBinning(nbins, xmin_xbj, xmax_xbj);
    std::vector<double> log_bins_x = CreateLogBinning(nbon, xmin_xbj, xM);
    std::vector<double> log_bins_Q = CreateLogBinning(nbon, qm, qM);
    TH2D *xQplane = new TH2D("xQplane", "Q^2 vs x_B | pion | y<0.95; x_B; Q^2", nbon, log_bins_x.data(), nbon, log_bins_Q.data());
    TH1D *truepionQ2 = new TH1D("truepionQ2", "Production of Pions+ in function of Q^2; GeV^2",  nbins, log_bins_Q2.data());
    TH1D *RecpionQ2 = new TH1D("RecpionQ2", "Reconstruction of Pions+ in function of Q^2; GeV^2",  nbins, log_bins_Q2.data());
    TH1D *truekaonQ2 = new TH1D("truekaonQ2", "Production of Kaons in function of Q^2; GeV^2",  nbins, log_bins_Q2.data());
    TH1D *ReckaonQ2 = new TH1D("ReckaonQ2", "Reconstruction of Kaons in function of Q^2; GeV^2",  nbins, log_bins_Q2.data());
    TH1D *trueprotonQ2 = new TH1D("trueprotonQ2", "Production of Protons in function of Q^2; GeV^2",  nbins, log_bins_Q2.data());
    TH1D *RecprotonQ2 = new TH1D("RecprotonQ2", "Production of Protons in function of Q^2; GeV^2",  nbins, log_bins_Q2.data());
    // GRAFICI DI x_B
    TH1D *truepion_xbj = new TH1D("truepion_xbj", "Production of Pions+ in function of x_Bj; x_Bj", nbins, log_bins_xbj.data());
    TH1D *Recpion_xbj = new TH1D("Recpion_xbj", "Reconstruction of Pions+ in function of x_Bj; x_Bj", nbins, log_bins_xbj.data());
    TH1D *truekaon_xbj = new TH1D("truekaon_xbj", "Production of Kaons in function of x_Bj;", nbins, log_bins_xbj.data());
    TH1D *Reckaon_xbj = new TH1D("Reckaon_xbj", "Reconstruction of Kaons in function of x_Bj;", nbins, log_bins_xbj.data());
    TH1D *trueproton_xbj = new TH1D("trueproton_xbj", "Production of Protons in function of x_Bj", nbins, log_bins_xbj.data());
    TH1D *Recproton_xbj = new TH1D("Recproton_xbj", "Reconstruction of Protons in function of x_Bj", nbins, log_bins_xbj.data());
    // GRAFICI DI z
    double xmin_z = 1e-4;
    double xmax_z = 1;
    std::vector<double> log_bins_z = CreateLogBinning(nbins, xmin_z, xmax_z);
    TH1D *truepion_z = new TH1D("truepion_z", "Production of Pions+ in function of z; z", nbins, log_bins_z.data());
    TH1D *Recpion_z = new TH1D("Recpion_z", "Reconstruction of Pions+ in function of z; z", nbins, log_bins_z.data());
    TH1D *truekaon_z = new TH1D("truekaon_z", "Production of Kaons in function of z", nbins, log_bins_z.data());
    TH1D *Reckaon_z = new TH1D("Reckaon_z", "Reconstruction of kaons in function of z", nbins, log_bins_z.data());
    TH1D *trueproton_z = new TH1D("trueproton_z", "Production of Protons in function of z", nbins, log_bins_z.data());
    TH1D *Recproton_z = new TH1D("Recproton_z", "Reconstruction of Protons in function of z", nbins, log_bins_z.data());
    // GRAFICI DI P_hT
    double xmin_PhT = 1e-2;
    double xmax_PhT = 50;
    std::vector<double> log_bins_PhT = CreateLogBinning(nbins, xmin_PhT, xmax_PhT);
    TH1D *truepion_PhT = new TH1D("truepion_PhT", "Production of Pions+ in function of P_hT; GeV", nbins, log_bins_PhT.data());
    TH1D *Recpion_PhT = new TH1D("Recpion_PhT", "Reconstruction of Pions+ in function of P_hT; GeV", nbins, log_bins_PhT.data());
    TH1D *truekaon_PhT = new TH1D("truekaon_PhT", "Production of Kaons in function of P_hT; GeV", nbins, log_bins_PhT.data());
    TH1D *Reckaon_PhT = new TH1D("Reckaon_PhT", "Reconstruction of Kaons in function of P_hT; GeV", nbins, log_bins_PhT.data());
    TH1D *trueproton_PhT = new TH1D("trueproton_PhT", "Production of Protons in function of P_hT; GeV", nbins, log_bins_PhT.data());
    TH1D *Recproton_PhT = new TH1D("Recproton_PhT", "Reconstruction of Protons in function of P_hT; GeV", nbins, log_bins_PhT.data());
    // GRACIFI DEL mom
    double xmin_mom = 1e-1;
    double xmax_mom = 50;
    std::vector<double> log_bins_mom = CreateLogBinning(nbins, xmin_mom, xmax_mom);
    TH1D *truepion_mom = new TH1D("truepion_mom", "Production of Pions+ in function of Mom; GeV", 80, log_bins_mom.data());
    TH1D *Recpion_mom = new TH1D("Recpion_mom", "Reconstruction of Pions+ in function of Mom; GeV",  80, log_bins_mom.data());
    TH1D *truekaon_mom = new TH1D("truekaon_mom", "Production of Kaons in function of Mom; GeV",  80, log_bins_mom.data());
    TH1D *Reckaon_mom = new TH1D("Reckaon_mom", "Reconstruction of Kaons in function of Mom; GeV",  80, log_bins_mom.data());
    TH1D *trueproton_mom = new TH1D("trueproton_mom", "Production of Protons in function of Mom; GeV",  80, log_bins_mom.data());
    TH1D *Recproton_mom = new TH1D("Recproton_mom", "Reconstruction of Protons in function of Mom; GeV",  80, log_bins_mom.data());
    // PROVA GRAFICI 2D
    double xmin_q2 = 1;
    double xmax_q2 = 100;
    std::vector<double> log_bins_Mom = CreateLogBinning(nben, xmin_mom, xmax_mom);
    std::vector<double> log_bins_Q22 = CreateLogBinning(nben, xmin_q2, xmax_q2);
    std::vector<double> log_bins_PhT2 = CreateLogBinning(nben, xmin_PhT, xmax_PhT);
    std::vector<double> log_bins_z2 = CreateLogBinning(nben, xmin_z, xmax_z);
    std::vector<double> log_bins_xbj2 = CreateLogBinning(nben, xmin_xbj, xmax_xbj);
    TH2D *pion_MomVsQ2 = new TH2D("pion_MomVsQ2", "Mom vs Q^2 | pion | y<0.95; Mom [GeV]; Q^2 [GeV^2]", nben, 0, 50, nben, log_bins_Q22.data());
    TH2D *pion_MomVsEta = new TH2D("pion_MomVsEta", "Mom vs Eta | pion | y<0.95; Eta; Mom [GeV]", nben, -3.5 ,3.5, nben, 0, 50);
    TH2D *pion_MomVsPhi = new TH2D("pion_MomVsTheta", "Mom vs Theta (Polar) | pion | y<0.95; Theta [Deg]; Mom [GeV]", nben, 0, 180, nben, 0, 50);
    TH2D *pion_MomVsTheta = new TH2D("pion_MomVsPhi", "Mom vs Phi (Azimuth) | pion | y<0.95; Phi [Deg]; Mom [GeV]", nben, -180, 180, nben, 0, 50);
    TH2D *pion_MomVsTheta_HP = new TH2D("pion_MomVsTheta_HP", "Mom vs Theta (Polar, Hadron plane) | pion | y<0.95; Theta [Deg]; Mom [GeV]", nben, 0, 180, nben, 0, 50);
    TH2D *pion_MomVsPhi_HP = new TH2D("pion_MomVsPhi_HP", "Mom vs Phi (Azimuth, Hadron plane) | pion | y<0.95; Phi [Deg]; Mom [GeV]", nben, -180, 180, nben, 0, 50);
    TH2D *pion_MomVsEta_rec = new TH2D("pion_MomVsEta_rec", "Mom vs Eta (MC Reco) | pion | y<0.95; Eta; Mom [GeV]", nben, -3.5 ,3.5, nben, 0, 50);
    TH2D *pion_MomVsPhi_rec = new TH2D("pion_MomVsTheta_rec", "Mom vs Theta (Polar) (MC Reco) | pion | y<0.95; Theta [Deg]; Mom [GeV]", nben, 0, 180, nben, 0, 50);
    TH2D *pion_MomVsTheta_rec = new TH2D("pion_MomVsPhi_rec", "Mom vs Phi (Azimuth) (MC Reco) | pion | y<0.95; Phi [Deg]; Mom [GeV]", nben, -180, 180, nben, 0, 50);
    TH2D *pion_PhTVsQ2 = new TH2D("pion_PhTVsQ2", "PhT vs Q^2 | pion | y<0.95; PhT [GeV]; Q^2 [GeV^2]", nben, 0, 20, nben, log_bins_Q22.data());
    TH2D *pion_PhTvsMom = new TH2D("pion_PhTvsMom", "PhT Vs Mom | pion | y<0.95; PhT [GeV]; Mom [GeV]",nben, 0, 20, nben, 0, 50);
    TH2D *pion_PhTvsZ = new TH2D("pion_PhTvsZ", "Z vs PhT | pion | y<0.95; z; PhT [GeV]", nben, 0, 1, nben, 0, 20);
    TH2D *pion_PhTvsEta = new TH2D("pion_PhTVsEta", "PhT vs Eta | pion | y<0.95; Eta; PhT [GeV]", nben, -3.5 ,3.5, nben, 0, 20);
    TH2D *pion_PhTvsPhi = new TH2D("pion_PhTvsTheta", "Theta (Polar) vs PhT | pion | y<0.95; Theta [Deg]; PhT [GeV]", nben, 0, 180, nben, 0, 20);
    TH2D *pion_PhTvsTheta = new TH2D("pion_PhTvsPhi", "Phi (Azimuth) vs PhT | pion | y<0.95; Phi [Deg]; PhT [GeV]", nben, -180, 180, nben, 0, 20);
    TH2D *pion_PhTvsTheta_HP = new TH2D("pion_PhTvsTheta_HP", "Theta (Polar, Hadron plane) vs PhT | pion | y<0.95; Phi [Deg]; PhT [GeV]", nben, 0, 180, nben, 0, 20);
    TH2D *pion_PhTvsPhi_HP = new TH2D("pion_PhTvsPhi_HP", "Phi (Azimuth, Hadron plane) vs PhT | pion | y<0.95; Theta [Deg]; PhT [GeV]", nben, -180, 180, nben, 0, 20);
    TH2D *pion_PhTvsMom_rec = new TH2D("pion_PhTvsMom_rec", "PhT Vs Mom (MC Reco) | pion | y<0.95; PhT [GeV]; Mom [GeV]",nben, 0, 20, nben, 0, 50);
    TH2D *pion_PhTvsZ_rec = new TH2D("pion_PhTvsZ_rec", "Z vs PhT (MC Reco) | pion | y<0.95; z; PhT [GeV]", nben, 0, 1, nben, 0, 20);
    TH2D *pion_PhTvsEta_rec = new TH2D("pion_PhTVsEta_rec", "PhT vs Eta (MC Reco) | pion | y<0.95; Eta; PhT [GeV]", nben, -3.5 ,3.5, nben, 0, 20);
    TH2D *pion_PhTvsPhi_rec = new TH2D("pion_PhTvsTheta_rec", "Theta (Polar) vs PhT (MC Reco) | pion | y<0.95; Theta [Deg]; PhT [GeV]", nben, 0, 180, nben, 0, 20);
    TH2D *pion_PhTvsTheta_rec = new TH2D("pion_PhTvsPhi_rec", "Phi (Azimuth) vs PhT (MC Reco) | pion | y<0.95; Phi [Deg]; PhT [GeV]", nben, -180, 180, nben, 0, 20);
    TH2D *pion_ThetaVsPhi = new TH2D("pion_PhiVsTheta", "Phi (Azimuth) vs Theta (Polar) | pion | y<0.95; Phi [Deg]; Theta [Deg]", nben, -180, 180, nben, 0, 180);
    TH2D *pion_Theta_HPVsPhi_HP = new TH2D("pion_Phi_HPVsTheta_HP", "Phi_HP (Azimuth) vs Theta_HP (Polar)  | pion | y<0.95; Phi [Deg]; Theta [Deg]", nben, -180, 180, nben, 0, 180);
    TH2D *pion_ThetaVsPhi_rec = new TH2D("pion_PhiVsTheta_rec", "Phi (Azimuth) vs Theta (Polar) (MC Reco) | pion | y<0.95; Phi [Deg]; Theta [deg]", nben, -180, 180, nben, 0, 180);
    TH2D *pion_ZvsQ2 = new TH2D("pion_ZvsQ2", "Z vs Q^2 | pion | y<0.95; z ; Q^2 [GeV^2]", nben, 0, 1, nben, log_bins_Q22.data());
    TH2D *pion_ZvsMom = new TH2D("pion_ZvsMom", "Z vs Mom | pion | y<0.95; z; Mom [GeV]", nben, 0, 1, nben, 0, 50);
    TH2D *pion_ZvsPhi = new TH2D("pion_ZVsTheta", "Z vs Theta (Polar) | pion | y<0.95; Theta [Deg]; z", nben, 0, 180, nben, 0, 1);
    TH2D *pion_ZvsTheta = new TH2D("pion_ZvsPhi", "Z vs Phi (Azimuth) | pion | y<0.95; Phi [Deg]; z", nben, -180, 180, nben, 0, 1);
    TH2D *pion_ZvsTheta_HP = new TH2D("pion_ZVsTheta_HP", "Z vs Theta (Polar, Hadron plane) | pion | y<0.95; Theta [Deg]; z", nben, 0, 180, nben, 0, 1);
    TH2D *pion_ZvsPhi_HP = new TH2D("pion_ZvsPhi_HP", "Z vs Phi (Azimuth, Hadron plane) | pion | y<0.95; Phi [Deg]; z", nben, -180, 180, nben, 0, 1);
    TH2D *pion_ZvsMom_rec = new TH2D("pion_ZvsMom_rec", "Z vs Mom (MC Reco) | pion | y<0.95; z; Mom [GeV]", nben, 0, 1, nben, 0, 50);
    TH2D *pion_ZvsPhi_rec = new TH2D("pion_ZVsTheta_rec", "Z vs Theta (Polar) (MC Reco) | pion | y<0.95; Theta [Deg]; z", nben, 0, 180, nben, 0, 1);
    TH2D *pion_ZvsTheta_rec = new TH2D("pion_ZvsPhi_rec", "Z vs Phi (Azimuth) (MC Reco) | pion | y<0.95; Phi [Deg]; z", nben, -180, 180, nben, 0, 1);
    TH2D *pion_XbVsMom = new TH2D("pion_XbVsMom", "Xb vs Mom | pion | y<0.95; xB; Mom [GeV]", nben, log_bins_xbj2.data(), nben, 0, 50);
    TH2D *pion_XbVsPhT = new TH2D("pion_XbVsPhT", "Xb vs PhT | pion | y<0.95; xB; PhT [GeV]", nben, log_bins_xbj2.data(), nben, 0, 20);
    TH2D *pion_XbVsZ = new TH2D("pion_XbVsZ", "Xb vs Z | pion | y<0.95; xB; z", nben, log_bins_xbj2.data(), nben, 0, 1);
    TH2D *pion_XbVsEta = new TH2D("pion_XbVsEta", "Xb vs Eta | pion | y<0.95; xB; Eta",nben, -3.5, 3.5, nben, log_bins_xbj2.data());
    TH2D *pion_XbVsPhi = new TH2D("pion_XbVsTheta", "Xb vs Theta (Polar) | pion | y<0.95; Theta [Deg]; xB", nben, 0, 180, nben, log_bins_xbj2.data());
    TH2D *pion_XbVsTheta = new TH2D("pion_XbvsPhi", "Xb vs Phi (Azimuth) | pion | y<0.95; Phi [Deg]; xB", nben, -180, 180, nben, log_bins_xbj2.data());
    TH2D *pion_XbVsTheta_HP = new TH2D("pion_XbVsTheta_HP", "Xb vs Theta (Polar, Hadron plane) | pion | y<0.95; Theta [Deg]; xB", nben, 0, 180, nben, log_bins_xbj2.data());
    TH2D *pion_XbVsPhi_HP = new TH2D("pion_XbvsPhi_HP", "Xb vs Phi (Azimuth, Hadron plane) | pion | y<0.95; Phi [Deg]; xB", nben, -180, 180, nben, log_bins_xbj2.data());
    TH2D *pion_XbVsMom_rec = new TH2D("pion_XbVsMom_rec", "Xb vs Mom (MC Reco) | pion | y<0.95; xB; Mom [GeV]", nben, log_bins_xbj2.data(), nben, 0, 50);
    TH2D *pion_XbVsPhT_rec = new TH2D("pion_XbVsPhT_rec", "Xb vs PhT (MC Reco) | pion | y<0.95; xB; PhT [GeV]", nben, log_bins_xbj2.data(), nben, 0, 20);
    TH2D *pion_XbVsZ_rec = new TH2D("pion_XbVsZ_rec", "Xb vs Z (MC Reco) | pion | y<0.95; xB; z", nben, log_bins_xbj2.data(), nben, 0, 1);
    TH2D *pion_XbVsPhi_rec = new TH2D("pion_XbVsTheta_rec", "Xb vs Theta (Polar) (MC Reco) | pion | y<0.95; Theta [Deg]; xB", nben, 0, 180, nben, log_bins_xbj2.data());
    TH2D *pion_XbVsTheta_rec = new TH2D("pion_XbVsPhi_rec", "Xb vs Phi (Azimuth) (MC Reco) | pion | y<0.95; Phi [Deg]; xB", nben, -180, 180, nben, log_bins_xbj2.data());
    TH2D *pion_XbVsQ2 = new TH2D("pion_XbVsQ2", "Xb vs Q^2 | pion | y<0.95; xB; Q^2 [GeV^2]", nben, log_bins_xbj2.data(), nben, log_bins_Q22.data());

    TH2D *kaon_MomVsEta = new TH2D("kaon_MomVsEta", "Mom vs Eta | kaon | y<0.95; Eta; Mom [GeV]", nben, -3.5 ,3.5, nben, 0, 50);
    TH2D *kaon_MomVsPhi = new TH2D("kaon_MomVsTheta", "Mom vs Theta (Polar) | kaon | y<0.95; Theta [Deg]; Mom [GeV]", nben, 0, 180, nben, 0, 50);
    TH2D *kaon_MomVsTheta = new TH2D("kaon_MomVsPhi", "Mom vs Phi (Azimuth) | kaon | y<0.95; Phi [Deg]; Mom [GeV]", nben, -180, 180, nben, 0, 50);
    TH2D *kaon_MomVsEta_rec = new TH2D("kaon_MomVsEta_rec", "Mom vs Eta (MC Reco) | kaon | y<0.95; Eta; Mom [GeV]", nben, -3.5 ,3.5, nben, 0, 50);
    TH2D *kaon_MomVsPhi_rec = new TH2D("kaon_MomVsTheta_rec", "Mom vs Theta (Polar) (MC Reco) | kaon | y<0.95; Theta [Deg]; Mom [GeV]", nben, 0, 180, nben, 0, 50);
    TH2D *kaon_MomVsTheta_rec = new TH2D("kaon_MomVsPhi_rec", "Mom vs Phi (Azimuth) (MC Reco) | kaon | y<0.95; Phi [Deg]; Mom [GeV]", nben, -180, 180, nben, 0, 50);
    TH2D *kaon_PhTvsMom = new TH2D("kaon_PhTvsMom", "PhT Vs Mom | kaon | y<0.95; PhT [GeV]; Mom [GeV]",nben, 0, 20, nben, 0, 50);
    TH2D *kaon_PhTvsZ = new TH2D("kaon_PhTvsZ", "Z vs PhT | kaon | y<0.95; z; PhT [GeV]", nben, 0, 1, nben, 0, 20);
    TH2D *kaon_PhTvsEta = new TH2D("kaon_PhTVsEta", "PhT vs Eta | kaon | y<0.95; Eta; PhT [GeV]", nben, -3.5 ,3.5, nben, 0, 20);
    TH2D *kaon_PhTvsPhi = new TH2D("kaon_PhTVsTheta", "Theta (Polar) vs PhT | kaon | y<0.95; Theta [Deg]; PhT [GeV]", nben, 0, 180, nben, 0, 20);
    TH2D *kaon_PhTvsTheta = new TH2D("kaon_PhTvsPhi", "Phi (Azimuth) vs PhT | kaon | y<0.95; Phi [Deg]; PhT [GeV]", nben, -180, 180, nben, 0, 20);
    TH2D *kaon_PhTvsMom_rec = new TH2D("kaon_PhTvsMom_rec", "PhT Vs Mom (MC Reco) | kaon | y<0.95; PhT [GeV]; Mom [GeV]",nben, 0, 20, nben, 0, 50);
    TH2D *kaon_PhTvsZ_rec = new TH2D("kaon_PhTvsZ_rec", "Z vs PhT (MC Reco) | kaon | y<0.95; z; PhT [GeV]", nben, 0, 1, nben, 0, 20);
    TH2D *kaon_PhTvsEta_rec = new TH2D("kaon_PhTVsEta_rec", "PhT vs Eta (MC Reco) | kaon | y<0.95; Eta; PhT [GeV]", nben, -3.5 ,3.5, nben, 0, 20);
    TH2D *kaon_PhTvsPhi_rec = new TH2D("kaon_PhTVsTheta_rec", "Theta (Polar) vs PhT (MC Reco) | kaon | y<0.95; Theta [Deg]; PhT [GeV]", nben, 0, 180, nben, 0, 20);
    TH2D *kaon_PhTvsTheta_rec = new TH2D("kaon_PhTvsPhi_rec", "Phi (Azimuth) vs PhT (MC Reco) | kaon | y<0.95; Phi [Deg]; PhT [GeV]", nben, -180, 180, nben, 0, 20);
    TH2D *kaon_ThetaVsPhi = new TH2D("kaon_PhiVsTheta", "Phi (Azimuth) vs Theta (Polar) | kaon | y<0.95; Phi [Deg]; Theta [Deg]", nben, -180, 180, nben, 0, 180);
    TH2D *kaon_ThetaVsPhi_rec = new TH2D("kaon_PhiVsTheta_rec", "Phi (Azimuth) vs Theta (Polar) (MC Reco) | kaon | y<0.95; Phi [Deg]; Theta [deg]", nben, -180, 180, nben, 0, 180);
    TH2D *kaon_ZvsMom = new TH2D("kaon_ZvsMom", "Z vs Mom | kaon | y<0.95; z; Mom [GeV]", nben, 0, 1, nben, 0, 50);
    TH2D *kaon_ZvsPhi = new TH2D("kaon_ZVsTheta", "Z vs Theta (Polar) | kaon | y<0.95; Theta [Deg]; z", nben, 0, 180, nben, 0, 1);
    TH2D *kaon_ZvsTheta = new TH2D("kaon_ZvsPhi", "Z vs Phi (Azimuth) | kaon | y<0.95; Phi [Deg]; z", nben, -180, 180, nben, 0, 1);
    TH2D *kaon_ZvsMom_rec = new TH2D("kaon_ZvsMom_rec", "Z vs Mom (MC Reco) | kaon | y<0.95; z; Mom [GeV]", nben, 0, 1, nben, 0, 50);
    TH2D *kaon_ZvsPhi_rec = new TH2D("kaon_ZVsTheta_rec", "Z vs Theta (Polar) (MC Reco) | kaon | y<0.95; Theta [Deg]; z", nben, 0, 180, nben, 0, 1);
    TH2D *kaon_ZvsTheta_rec = new TH2D("kaon_ZvsPhi_rec", "Z vs Phi (Azimuth) (MC Reco) | kaon | y<0.95; Phi [Deg]; z", nben, -180, 180, nben, 0, 1);
    TH2D *kaon_XbVsMom = new TH2D("kaon_XbVsMom", "Xb vs Mom | kaon | y<0.95; xB; Mom [GeV]", nben, log_bins_xbj2.data(), nben, 0, 50);
    TH2D *kaon_XbVsPhT = new TH2D("kaon_XbVsPhT", "Xb vs PhT | kaon | y<0.95; xB; PhT [GeV]", nben, log_bins_xbj2.data(), nben, 0, 20);
    TH2D *kaon_XbVsZ = new TH2D("kaon_XbVsZ", "Xb vs Z | kaon | y<0.95; xB; z", nben, log_bins_xbj2.data(), nben, 0, 1);
    TH2D *kaon_XbVsPhi = new TH2D("kaon_XbVsTheta", "Xb vs Theta (Polar) | kaon | y<0.95; Theta [Deg]; xB", nben, 0, 180, nben, log_bins_xbj2.data());
    TH2D *kaon_XbVsTheta = new TH2D("kaon_XbvsPhi", "Xb vs Phi (Azimuth) | kaon | y<0.95; Phi [Deg]; xB", nben, -180, 180, nben, log_bins_xbj2.data());
    TH2D *kaon_XbVsMom_rec = new TH2D("kaon_XbVsMom_rec", "Xb vs Mom (MC Reco) | kaon | y<0.95; xB; Mom [GeV]", nben, log_bins_xbj2.data(), nben, 0, 50);
    TH2D *kaon_XbVsPhT_rec = new TH2D("kaon_XbVsPhT_rec", "Xb vs PhT (MC Reco) | kaon | y<0.95; xB; PhT [GeV]", nben, log_bins_xbj2.data(), nben, 0, 20);
    TH2D *kaon_XbVsZ_rec = new TH2D("kaon_XbVsZ_rec", "Xb vs Z (MC Reco) | kaon | y<0.95; xB; z", nben, log_bins_xbj2.data(), nben, 0, 1);
    TH2D *kaon_XbVsPhi_rec = new TH2D("kaon_XbVsTheta_rec", "Xb vs Theta (Polar) (MC Reco) | kaon | y<0.95; Theta [Deg]; xB", nben, 0, 180, nben, log_bins_xbj2.data());
    TH2D *kaon_XbVsTheta_rec = new TH2D("kaon_XbVsPhi_rec", "Xb vs Phi (Azimuth) (MC Reco) | kaon | y<0.95; Phi [Deg]; xB", nben, -180, 180, nben, log_bins_xbj2.data());

    TH2D *proton_MomVsEta = new TH2D("proton_MomVsEta", "Mom vs Eta | proton | y<0.95; Eta; Mom [GeV]", nben, -3.5 ,3.5, nben, 0, 50);
    TH2D *proton_MomVsPhi = new TH2D("proton_MomVsTheta", "Mom vs Theta (Polar) | proton | y<0.95; Theta [Deg]; Mom [GeV]", nben, 0, 180, nben, 0, 50);
    TH2D *proton_MomVsEta_rec = new TH2D("proton_MomVsEta_rec", "Mom vs Eta | proton | y<0.95; Eta; Mom [GeV]", nben, -3.5 ,3.5, nben, 0, 50);
    TH2D *proton_MomVsPhi_rec = new TH2D("proton_MomVsTheta_rec", "Mom vs Theta (Polar) | proton | y<0.95; Theta [Deg]; Mom [GeV]", nben, 0, 180, nben, 0, 50);
    
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
                        20, 0, 180, 20, -180, 180);
    TH2D *h3 = new TH2D("h3", "Proiezione sferica; Theta (deg); Phi (deg); Intensità", 
                        20, 0, 180, 20, -180, 180);
    //PROVA DI TH3D
    int nban = 10;
    std::vector<double> log_bins_Mom3 = CreateLogBinning(nban, xmin_mom, xmax_mom);
    std::vector<double> log_bins_PhT3 = CreateLogBinning(nban, xmin_PhT, xmax_PhT);
    std::vector<double> log_bins_z3 = CreateLogBinning(nban, xmin_z, xmax_z);
    std::vector<double> log_bins_xbj3 = CreateLogBinning(nban, xmin_xbj, xmax_xbj);
    TH3D *pion_MomVsPhiVsTheta = new TH3D("pion_MomVsPhiVsTheta", 
    "Mom vs Theta (Polar) vs Phi (Azimuth) | pion | y<0.95; Theta [Deg]; Mom [GeV]; Phi [Deg]", 
    nban, -180, 180, nban, 0, 1, nban, 0, 180);
    pion_MomVsPhiVsTheta->GetYaxis()->Set(nban, log_bins_Mom3.data());
    TH3D *pion_PhTvsZvsXb = new TH3D("pion_PhTvsZvsXb",
    "PhT vs Z vs xB | pion | y<0.95; xB; z; PhT [GeV]",
    nban, 0, 1, nban, 0, 1, nban, 0, 100);
    pion_PhTvsZvsXb->GetXaxis()->Set(nban, log_bins_xbj3.data());
    pion_PhTvsZvsXb->GetYaxis()->Set(nban, log_bins_z3.data());
    pion_PhTvsZvsXb->GetZaxis()->Set(nban, log_bins_PhT3.data());
    TH3D *z_axisSP = new TH3D("z_axisSP", "Scattered photon polar angle; x; y; z",  
    30, 0, 1, 30, 0, 1, 30, 0, 1);

      
    // "REAL"
    TH1D *RealpionQ2 = new TH1D("RealpionQ2", "Real reconstruction of Pions+ in function of Q^2; GeV^2",  nbins, log_bins_Q2.data());
    TH1D *RealkaonQ2 = new TH1D("RealkaonQ2", "Reconstruction of Kaons in function of Q^2; GeV^2",  nbins, log_bins_Q2.data());
    TH1D *RealprotonQ2 = new TH1D("RealprotonQ2", "Production of Protons in function of Q^2; GeV^2",  nbins, log_bins_Q2.data());
    TH1D *RealpionX = new TH1D("RealpionX", "Real reconstruction of Pions+ in function of Q^2; GeV^2",  nbins, log_bins_xbj.data());
    TH1D *RealkaonX = new TH1D("RealkaonX", "Reconstruction of Kaons in function of Q^2; GeV^2",  nbins, log_bins_xbj.data());
    TH1D *RealprotonX = new TH1D("RealprotonX", "Production of Protons in function of Q^2; GeV^2",  nbins, log_bins_xbj.data());
    TH1D *RealpionZ = new TH1D("RealpionZ", "Real reconstruction of Pions+ in function of Q^2; GeV^2",  nbins, log_bins_z.data());
    TH1D *RealkaonZ = new TH1D("RealkaonZ", "Reconstruction of Kaons in function of Q^2; GeV^2",  nbins, log_bins_z.data());
    TH1D *RealprotonZ = new TH1D("RealprotonZ", "Production of Protons in function of Q^2; GeV^2",  nbins, log_bins_z.data());
    TH1D *RealpionPhT = new TH1D("RealpionPht", "Real reconstruction of Pions+ in function of Q^2; GeV^2",  nbins, log_bins_PhT.data());
    TH1D *RealkaonPhT= new TH1D("RealkaonPhT", "Reconstruction of Kaons in function of Q^2; GeV^2",  nbins, log_bins_PhT.data());
    TH1D *RealprotonPhT = new TH1D("RealprotonPhT", "Production of Protons in function of Q^2; GeV^2",  nbins, log_bins_PhT.data());
    TH1D *RealpionEta = new TH1D("RealpionEta", "Real reconstruction of Pions+ in function of Q^2; GeV^2",  80, -3.5 ,3.5);
    TH1D *RealkaonEta = new TH1D("RealkaonEta", "Reconstruction of Kaons in function of Q^2; GeV^2",  80, -3.5 ,3.5);
    TH1D *RealprotonEta = new TH1D("RealprotonEta", "Production of Protons in function of Q^2; GeV^2",  80, -3.5 ,3.5);
    TH1D *RealpionPhi = new TH1D("RealpionPhi", "Real reconstruction of Pions+ in function of Q^2; GeV^2",  80, 0, 180);
    TH1D *RealkaonPhi = new TH1D("RealkaonPhi", "Reconstruction of Kaons in function of Q^2; GeV^2",  80, 0, 180);
    TH1D *RealprotonPhi = new TH1D("RealprotonPhi", "Production of Protons in function of Q^2; GeV^2",  80, 0, 180);
    TH1D *RealpionMom = new TH1D("RealpionMom", "Real reconstruction of Pions+ in function of Q^2; GeV^2",  nbins, log_bins_mom.data());
    TH1D *RealkaonMom = new TH1D("RealkaonMom", "Reconstruction of Kaons in function of Q^2; GeV^2",  nbins, log_bins_mom.data());
    TH1D *RealprotonMom = new TH1D("RealprotonMom", "Production of Protons in function of Q^2; GeV^2",  nbins, log_bins_mom.data());

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
    TH1D *GammaPolarDeg = new TH1D("GammaPolarDeg", "Scattered photon polar angle; Theta [Deg]",  nben, 0, 180);
    TH1D *GammaAzimuthDeg = new TH1D("GammaAzimuthDeg", "Scattered photon azimuthal angle; Phi [Deg]",  nben, -180, 180);


    // ALCUNI VETTORI UTILI
    std::vector<TLorentzVector> scatElectron;
    std::vector<TVector3> recScatElectron;
    std::vector<TVector3> GammaVector;
    std::vector<TVector3> BeamElectronVector;
    //std::vector<TVector3> z_axis_SP 
    //std::vector<TVector3> ipsilon;
    //std::vector<TVector3> y_axis_SP;
    //std::vector<TVector3> x_axis_SP;
    std::vector<float> scatPhi;
    std::vector<float> PolarGammaRad;
    std::vector<float> PolarGammaDeg;
    std::vector<float> AzimuthGammaRad;
    std::vector<float> AzimuthGammaDeg;
    std::vector<TVector3> elMom_pion;
    TVector3 ElBeam(0.,0.,-18.); 
    TLorentzVector ElectronBeam(18., 0, 0, -18.);
    TLorentzVector ProtonBeam(275., 0, 0, 275.);
    std::vector<TLorentzVector> q;
    double currentPhi = 0;
    double currentMom = 0;
    std::vector<float> scatElPhipion;
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

    while(tree_reader.Next()) { // Loop over events

      for(unsigned int i=0; i<partGenStat.GetSize(); i++) // Loop over thrown particles
        {
          count += 1;
          TVector3 part(partMomX[i],partMomY[i],partMomZ[i]);
          float partEta = part.PseudoRapidity();
          double phis = part.Theta();
          double y_DA= ((TMath::Tan(phis*0.5)) / (TMath::Tan(phis*0.5) + TMath::Tan(currentPhi*0.5)));
          if(partEta <= 3.5){
            if(partEta >= -3.5){
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
                float Pmom = BeamProton.Mag();
                float Peng = sqrt(Pmom*Pmom + 0.939*0.939);
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
              }
              if(pdg == 11){
                TVector3 BeamElectron(partMomX[i],partMomY[i],partMomZ[i]);
                float Emom = BeamElectron.Mag();
                double Beam_EThetaRad = BeamElectron.Phi();
                double Beam_EThetaDeg = Beam_EThetaRad * (180.0/TMath::Pi());
                double Beam_EPhiRad = BeamElectron.Theta();
                double Beam_EPhiDeg = Beam_EPhiRad * (180.0/TMath::Pi());
                double Ex = std::abs(1) * TMath::Sin(Beam_EPhiRad) * TMath::Cos(Beam_EThetaRad);
                double Ey = std::abs(1) * TMath::Sin(Beam_EPhiRad) * TMath::Sin(Beam_EThetaRad);
                double Ez = std::abs(1) * TMath::Cos(Beam_EPhiRad);

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
                  float mom = ElMom.Mag();
                  double etta = ElMom.PseudoRapidity();
                  if(std::abs(etta) <= 3.5){
                    if(parentsIndex[i]<=500){
                      TLorentzVector tlv(mom, ElMom.X(), ElMom.Y(), ElMom.Z());
                      float angleR = ElMom.Theta();
                      float angle = angleR * (180.0 / TMath::Pi());
                      currentPhi = ElMom.Theta();
                      currentMom = ElMom.Mag();
                      currentQ2pion.SetXYZ(partMomX[i],partMomY[i],partMomZ[i]);
                      currentQ2kaon.SetXYZ(partMomX[i],partMomY[i],partMomZ[i]);
                      currentQ2proton.SetXYZ(partMomX[i],partMomY[i],partMomZ[i]);

                      scatEl->Fill(mom);
                      scatElAngle->Fill(angle);
                      //scatElectron.push_back(ElMom);  // to use it outside the cycle
                      scatElectron.push_back(tlv);
                      scatPhi.push_back(angle);
                      BeamElectronVector.push_back(ElMom);

                        for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
                        {
                          if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
                            {
                              TVector3 recElmom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); 
                              float momE = recElmom.Mag();
                              recScatEl->Fill(momE);
                              recScatElectron.push_back(recElmom);
                              
                            }
                        }
                    }
                  }
                }
                if(pdg == 22){ // stiamo cercando il fotone scatterato
                  if(parentsIndex[i]<=6){
                    TVector3 ScatPhoton(partMomX[i],partMomY[i],partMomZ[i]);
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
                    GammaVector.push_back(ScatPhoton);
                  }
                }
            }
          // IDK IF THIS WILL WORK HERE...
          for (const auto& vec : scatElectron) {
            q.push_back(ElectronBeam - vec);
          }
          if(partGenStat[i] == 1) // Select stable thrown particles
            {
              //int pdg = partPdg[i];
              //if(pdg == 11 || pdg == 13 || pdg == 211 || pdg == 321 || pdg == 2212) // Look at charged particles (electrons, muons, pions, kaons, protons)  
              // PION   
              if(pdg == 211)
                { 
                  TVector3 truePionMom(partMomX[i],partMomY[i],partMomZ[i]);
                  float mom_pion = truePionMom.Mag();
                  float E_pion = sqrt(mom_pion*mom_pion + 0.139*0.139);
                  TLorentzVector pion(E_pion, partMomX[i],partMomY[i],partMomZ[i]);
                  float pionEta = truePionMom.PseudoRapidity();

                  if(pionEta >= -3.5 && pionEta <= 3.5){
                    double pionPhiRad = truePionMom.Theta(); // MOLTO CONFUSO LO SO, PHI (MIO) ANGOLO POLARE
                    double pionThetaRad = truePionMom.Phi() + TMath::Pi(); // THETA ANGOLO AZIMUTALE 
                    if(pionThetaRad >= TMath::Pi()){
                      pionThetaRad -= 2*TMath::Pi();
                    }
                    float pionPhi = pionPhiRad * (180.0 / TMath::Pi());
                    float pionTheta = pionThetaRad * (180.0 / TMath::Pi());
                    TLorentzVector scatElpion(currentMom, currentQ2pion.X(), currentQ2pion.Y(), currentQ2pion.Z());
                    TLorentzVector photon_pion = ElectronBeam - scatElpion;                

                    for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
                      {
                        int recpdg = std::abs(recPdg[j]);
                        if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
                          {
                            TVector3 recPionMom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); 
                            float mom_pion_Rec = recPionMom.Mag();
                            float E_pion_rec = sqrt(mom_pion_Rec*mom_pion_Rec + 0.139*0.139);
                            TLorentzVector pion_Rec(E_pion_rec, trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
                            float recpionEta = recPionMom.PseudoRapidity();

                            if(recpionEta <= 3.5 && recpionEta >= -3.5){
                              double recpPhi = recPionMom.Theta();
                              double recpTheta = recPionMom.Phi() + TMath::Pi();
                              if(recpTheta >= TMath::Pi()){
                                recpTheta -= 2*TMath::Pi();
                              }
                              float recpionTheta = recpTheta * (180.0 / TMath::Pi());
                              float recpionPhi = recpPhi * (180.0 / TMath::Pi());
                              TLorentzVector scatElpion_Rec(currentMom, currentQ2pion.X(), currentQ2pion.Y(), currentQ2pion.Z());
                              TLorentzVector photon_pion_Rec = ElectronBeam - scatElpion_Rec;

                              double y_DA_pion_Rec = ((TMath::Tan(recpPhi*0.5)) / (TMath::Tan(recpPhi*0.5) + TMath::Tan(currentPhi*0.5)));
                              if(y_DA_pion_Rec <= 0.95){
                                if(parentsIndex[i]<=5){
                                  double Q2_DA_pion_Rec = ((4*18*18*(1-y_DA_pion_Rec)) / (TMath::Tan(currentPhi*0.5)*TMath::Tan(currentPhi*0.5)));
                                  double xbj_DA_pion_Rec = (Q2_DA_pion_Rec) / (4*18*275*y_DA_pion_Rec);
                                  double z_DA_pion_Rec = (ProtonBeam * pion_Rec) / (ProtonBeam * photon_pion_Rec);
                                  TVector3 PhT_pion_vec_Rec = recPionMom - ((recPionMom * currentQ2pion) / currentQ2pion.Mag())*currentQ2pion.Unit();
                                  double PhT_pion_Rec = PhT_pion_vec_Rec.Mag();
                                  if(z_DA_pion_Rec <= 1){
                                    RecpionQ2->Fill(Q2_DA_pion_Rec);
                                    Recpion_xbj->Fill(xbj_DA_pion_Rec);
                                    Recpion_z->Fill(z_DA_pion_Rec);
                                    Recpion_PhT->Fill(PhT_pion_Rec);
                                    Recpion_mom->Fill(mom_pion_Rec);
                                    RecpionEta->Fill(recpionEta);
                                    RecpionPhi->Fill(recpionPhi);
                                    pion_MomVsEta_rec->Fill(recpionEta, mom_pion_Rec);
                                    pion_PhTvsEta_rec->Fill(recpionEta, PhT_pion_Rec);
                                    pion_MomVsPhi_rec->Fill(recpionPhi, mom_pion_Rec);
                                    pion_PhTvsZ_rec->Fill(z_DA_pion_Rec, PhT_pion_Rec);
                                    pion_MomVsTheta_rec->Fill(recpionTheta, mom_pion_Rec);
                                    pion_ThetaVsPhi_rec->Fill(recpionTheta,recpionPhi);
                                    pion_PhTvsPhi_rec->Fill(recpionPhi, PhT_pion_Rec);
                                    pion_PhTvsTheta_rec->Fill(recpionTheta, PhT_pion_Rec);
                                    pion_PhTvsMom_rec->Fill(PhT_pion_Rec, mom_pion_Rec);
                                    pion_ZvsMom_rec->Fill(z_DA_pion_Rec, mom_pion_Rec);
                                    pion_ZvsPhi_rec->Fill(recpionPhi,z_DA_pion_Rec);
                                    pion_ZvsTheta_rec->Fill(recpionTheta, z_DA_pion_Rec);
                                    pion_XbVsMom_rec->Fill(xbj_DA_pion_Rec, mom_pion_Rec);
                                    pion_XbVsPhT_rec->Fill(xbj_DA_pion_Rec, PhT_pion_Rec);
                                    pion_XbVsZ_rec->Fill(xbj_DA_pion_Rec, z_DA_pion_Rec);
                                    pion_XbVsPhi_rec->Fill(recpionPhi, xbj_DA_pion_Rec);
                                    pion_XbVsTheta_rec->Fill(recpionTheta, xbj_DA_pion_Rec);

                                    if(recpdg==pdg){
                                      RealpionQ2->Fill(Q2_DA_pion_Rec);
                                      RealpionX->Fill(xbj_DA_pion_Rec);
                                      RealpionZ->Fill(z_DA_pion_Rec);
                                      RealpionPhT->Fill(PhT_pion_Rec);
                                      RealpionEta->Fill(recpionEta);
                                      RealpionPhi->Fill(recpionPhi);
                                      RealpionMom->Fill(mom_pion_Rec);
                                    }
                                  }
                                  PIDgen2rec->Fill(kpion, checkRecoID(recpdg));
                                }
                              }
                            }
                          }
                      }
                    scatElPhipion.push_back(currentPhi); // that generate an array with the angle of the scat. el with the # of pions
                    scatElq_pion.push_back(currentQ2pion);
                    double y_DA_pion = ((TMath::Tan(pionPhiRad*0.5)) / (TMath::Tan(pionPhiRad*0.5) + TMath::Tan(currentPhi*0.5)));
                    // PROVIAMO A METTERE LE COSE DEL FOTONE
                    // prendo le tre componenti del fotone, siccome i loro versori sono allineati con gli assi del piano di scattering
                    TVector3 CurrentGammaVector = GammaVector[i];
                    TVector3 z_axis_SP = CurrentGammaVector.Unit(); // voglio che l'asse z sia lungo la direzione di gamma
                    double z_axis_z = z_axis_SP.Z();
                    double z_axis_y = z_axis_SP.Y();
                    double z_axis_x = z_axis_SP.X();
                    // per calcolare y mi serve il prodotto vettoriale tra il leptone entrante e gamma
                    TVector3 CurrentBeamElectronVector = BeamElectronVector[i];
                    TVector3 ipsilon = CurrentGammaVector.Cross(CurrentBeamElectronVector);
                    TVector3 y_axis_SP = ipsilon.Unit();
                    TVector3 x_axis_SP = y_axis_SP.Cross(z_axis_SP);
                    // BENE ADESSO ABBIAMO I NOSTRI ASSI DEL PIANO DI SCATTERING
                    double P_T = truePionMom.Perp(z_axis_SP); // questo comando mi da il vettore ortogonale a z, quindi il momento trasverso
                    double Theta_HP_Rad = truePionMom.Angle(z_axis_SP); // dovrebbe fornire l'angolo polare del hadron plane
                    double Theta_HP_Deg = Theta_HP_Rad * (180.0/TMath::Pi());
                    // però vogliamo anche il vettore del momento trasverso
                    double PdotZ = truePionMom.Dot(z_axis_SP); // momento proiettato sull'asse z
                    TVector3 Momentum_Z = PdotZ * z_axis_SP;
                    TVector3 P_T_Vector = truePionMom - Momentum_Z; // ora lo possiamo usare per trovare l'angolo azimutale
                    double P_T_Vector_x = P_T_Vector * x_axis_SP;
                    double P_T_Vector_y = P_T_Vector * y_axis_SP; // calcolo il momento trasverso sugli assi x e y
                    double Phi_HP_Rad = std::atan2(P_T_Vector_y, P_T_Vector_x);
                    double Phi_HP_Deg = Phi_HP_Rad * (180.0/TMath::Pi()); 

                    if(y_DA_pion <= 0.95){  
                        double Q2_DA_pion = ((4*18*18*(1-y_DA_pion)) / (TMath::Tan(currentPhi*0.5)*TMath::Tan(currentPhi*0.5)));
                        //double Q2_DA_pion = (scp*scp*TMath::Sin(currentPhi)*TMath::Sin(currentPhi))/(1-y_DA_pion);
                        double xbj_DA_pion = (Q2_DA_pion) / (4*18*275*y_DA_pion);
                        double z_DA_pion = (ProtonBeam * pion) / (ProtonBeam * photon_pion);
                        TVector3 PhT_pion_vec = truePionMom - ((truePionMom * currentQ2pion) / currentQ2pion.Mag())*currentQ2pion.Unit();
                        double PhT_pion = PhT_pion_vec.Mag();
                        TLorentzVector sq = ElectronBeam + ProtonBeam;
                        if(parentsIndex[i]<=5 ){
                            if(z_DA_pion <= 1){
                                truepionEta->Fill(pionEta);
                                truepionPhi->Fill(pionPhi);
                                truepionQ2->Fill(Q2_DA_pion); 
                                truepion_xbj->Fill(xbj_DA_pion);
                                truepion_z->Fill(z_DA_pion);
                                truepion_PhT->Fill(PhT_pion);
                                truepion_mom->Fill(mom_pion);
                                xQplane->Fill(xbj_DA_pion, Q2_DA_pion);
                                pion_MomVsQ2->Fill(mom_pion, Q2_DA_pion);
                                pion_MomVsEta->Fill(pionEta, mom_pion);
                                pion_PhTVsQ2->Fill(PhT_pion, Q2_DA_pion);
                                pion_PhTvsEta->Fill(pionEta, PhT_pion);
                                pion_MomVsPhi->Fill(pionPhi, mom_pion);
                                pion_PhTvsZ->Fill(z_DA_pion, PhT_pion);
                                pion_MomVsTheta->Fill(pionTheta, mom_pion);
                                pion_MomVsTheta_HP->Fill(Theta_HP_Deg, mom_pion);
                                pion_MomVsPhi_HP->Fill(Phi_HP_Deg, mom_pion);
                                pion_ThetaVsPhi->Fill(pionTheta,pionPhi);
                                pion_Theta_HPVsPhi_HP->Fill(Phi_HP_Deg, Theta_HP_Deg);
                                pion_PhTvsPhi->Fill(pionPhi, PhT_pion);
                                pion_PhTvsTheta->Fill(pionTheta, PhT_pion);
                                pion_PhTvsPhi_HP->Fill(Phi_HP_Deg, PhT_pion);
                                pion_PhTvsTheta_HP->Fill(Theta_HP_Deg, PhT_pion);
                                pion_MomVsPhiVsTheta->Fill(pionTheta, mom_pion, pionPhi);
                                pion_PhTvsMom->Fill(PhT_pion,mom_pion);
                                pion_ZvsQ2->Fill(z_DA_pion, Q2_DA_pion);
                                pion_ZvsMom->Fill(z_DA_pion, mom_pion);
                                pion_ZvsPhi->Fill(pionPhi,z_DA_pion);
                                pion_ZvsTheta->Fill(pionTheta, z_DA_pion);
                                pion_ZvsPhi_HP->Fill(Phi_HP_Deg,z_DA_pion);
                                pion_ZvsTheta_HP->Fill(Theta_HP_Deg, z_DA_pion);
                                pion_XbVsMom->Fill(xbj_DA_pion, mom_pion);
                                pion_XbVsPhT->Fill(xbj_DA_pion, PhT_pion);
                                pion_XbVsZ->Fill(xbj_DA_pion, z_DA_pion);
                                pion_XbVsPhi->Fill(pionPhi, xbj_DA_pion);
                                pion_XbVsTheta->Fill(pionTheta, xbj_DA_pion);
                                pion_XbVsPhi_HP->Fill(Phi_HP_Deg, xbj_DA_pion);
                                pion_XbVsTheta_HP->Fill(Theta_HP_Deg, xbj_DA_pion);
                                pion_XbVsQ2->Fill(xbj_DA_pion, Q2_DA_pion);
                                pion_PhTvsZvsXb->Fill(xbj_DA_pion, z_DA_pion, PhT_pion);
                                pion_XbVsEta->Fill(pionEta, xbj_DA_pion);
                                z_axisSP->Fill(z_axis_x, z_axis_y, z_axis_z);
                                /*
                                if(xbj_DA_pion<=1){
                                  for(int ix = 1; ix <= nben; ++ix){
                                    double uVal = pdf->xfxQ(2, xbj_DA_pion, Q2_DA_pion); 
                                    hPDF->SetBinContent(i, ix, uVal);
                                  }
                                  if(Q2_DA_pion >= 3 && Q2_DA_pion <= 3.2){
                                    pion_count += 1;
                                    double uVal2 = pdf->xfxQ(2, xbj_DA_pion, Q2_DA_pion);
                                    double dVal2 = pdf->xfxQ(1, xbj_DA_pion, Q2_DA_pion);
                                    double glVal2 = pdf->xfxQ(21, xbj_DA_pion, Q2_DA_pion);
                                    val_x.push_back(xbj_DA_pion);
                                    if(uVal2 >0 && uVal2 <= 1){
                                      uPDF->SetPoint(i, xbj_DA_pion, uVal2);
                                      val_uPDF.push_back(uVal2);
                                    }
                                    if(dVal2 > 0 && dVal2 < 1){
                                      dPDF->SetPoint(i, xbj_DA_pion, dVal2);
                                      val_dPDF.push_back(dVal2);
                                    }
                                    if(glVal2 > 0 && glVal2 >= 1){
                                      glPDF->SetPoint(i, xbj_DA_pion, 0.1*glVal2);
                                      val_glPDF.push_back(glVal2);
                                    }
                                      
                                  }
                                } */
                            }
                        }
                    }
                  }
                }  
              // KAON
              if(pdg == 321)
                {
                  TVector3 trueKaonMom(partMomX[i],partMomY[i],partMomZ[i]);
                  float mom_kaon = trueKaonMom.Mag();
                  float E_kaon = sqrt(mom_kaon*mom_kaon + 0.493*0.493);
                  TLorentzVector kaon(E_kaon, partMomX[i],partMomY[i],partMomZ[i]);
                  float kaonEta = trueKaonMom.PseudoRapidity();

                  if(kaonEta <= 3.5 && kaonEta >= -3.5){
                    double kaonPhiRad = trueKaonMom.Theta();
                    double kaonThetaRad = trueKaonMom.Phi() + TMath::Pi(); // THETA ANGOLO AZIMUTALE 
                    if(kaonThetaRad >= TMath::Pi()){
                      kaonThetaRad -= 2*TMath::Pi();
                    }
                    float kaonPhi = kaonPhiRad * (180.0 / TMath::Pi());
                    float kaonTheta = kaonThetaRad * (180.0 / TMath::Pi());
                    TLorentzVector scatElkaon(currentMom, currentQ2kaon.X(), currentQ2kaon.Y(), currentQ2kaon.Z());
                    TLorentzVector photon_kaon = ElectronBeam - scatElkaon;

                    for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
                      {
                        int recpdg = std::abs(recPdg[j]);
                        if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
                          {
                            TVector3 recKaonMom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); 
                            float mom_kaon_Rec = recKaonMom.Mag();
                            float E_kaon_rec = sqrt(mom_kaon_Rec*mom_kaon_Rec);
                            TLorentzVector kaon_Rec(E_kaon_rec, trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
                            float reckaonEta = recKaonMom.PseudoRapidity();

                            if(reckaonEta <= 3.5 && reckaonEta >= -3.5){
                              double reckPhi = recKaonMom.Theta();
                              double reckTheta = recKaonMom.Phi() + TMath::Pi();
                              if(reckTheta >= TMath::Pi()){
                                reckTheta -= 2*TMath::Pi();
                              }
                              float reckaonTheta = reckTheta * (180.0 / TMath::Pi());
                              float reckaonPhi = reckPhi * (180.0 / TMath::Pi());
                              TLorentzVector scatElkaon_Rec(currentMom, currentQ2kaon.X(), currentQ2kaon.Y(), currentQ2kaon.Z());
                              TLorentzVector photon_kaon_Rec = ElectronBeam - scatElkaon_Rec;
                              double y_DA_kaon_Rec = ((TMath::Tan(reckPhi*0.5)) / (TMath::Tan(reckPhi*0.5) + TMath::Tan(currentPhi*0.5)));
                              if(y_DA_kaon_Rec <= 0.95){
                                if(parentsIndex[i]<=5){
                                  double Q2_DA_kaon_Rec = ((4*18*18*(1-y_DA_kaon_Rec)) / (TMath::Tan(currentPhi*0.5)*TMath::Tan(currentPhi*0.5)));
                                  double xbj_DA_kaon_Rec = (Q2_DA_kaon_Rec) / (4*18*275*y_DA_kaon_Rec);
                                  double z_DA_kaon_Rec = (ProtonBeam * kaon_Rec) / (ProtonBeam * photon_kaon_Rec);
                                  TVector3 PhT_kaon_vec_Rec = recKaonMom - ((recKaonMom * currentQ2kaon) / currentQ2kaon.Mag())*currentQ2kaon.Unit();
                                  double PhT_kaon_Rec = PhT_kaon_vec_Rec.Mag();
                                  if(z_DA_kaon_Rec <= 1){
                                    ReckaonQ2->Fill(Q2_DA_kaon_Rec);
                                    Reckaon_xbj->Fill(xbj_DA_kaon_Rec);
                                    Reckaon_z->Fill(z_DA_kaon_Rec);
                                    Reckaon_PhT->Fill(PhT_kaon_Rec);
                                    Reckaon_mom->Fill(mom_kaon_Rec);
                                    ReckaonEta->Fill(reckaonEta);
                                    ReckaonPhi->Fill(reckaonPhi);
                                    kaon_MomVsEta_rec->Fill(reckaonEta, mom_kaon_Rec);
                                    kaon_PhTvsEta_rec->Fill(reckaonEta, PhT_kaon_Rec);
                                    kaon_MomVsPhi_rec->Fill(reckaonPhi, mom_kaon_Rec);
                                    kaon_PhTvsZ_rec->Fill(z_DA_kaon_Rec, PhT_kaon_Rec);
                                    kaon_MomVsTheta_rec->Fill(reckaonTheta, mom_kaon_Rec);
                                    kaon_ThetaVsPhi_rec->Fill(reckaonTheta,reckaonPhi);
                                    kaon_PhTvsPhi_rec->Fill(reckaonPhi, PhT_kaon_Rec);
                                    kaon_PhTvsTheta_rec->Fill(reckaonTheta, PhT_kaon_Rec);
                                    kaon_PhTvsMom_rec->Fill(PhT_kaon_Rec, mom_kaon_Rec);
                                    kaon_ZvsMom_rec->Fill(z_DA_kaon_Rec, mom_kaon_Rec);
                                    kaon_ZvsPhi_rec->Fill(reckaonPhi,z_DA_kaon_Rec);
                                    kaon_ZvsTheta_rec->Fill(reckaonTheta, z_DA_kaon_Rec);
                                    kaon_XbVsMom_rec->Fill(xbj_DA_kaon_Rec, mom_kaon_Rec);
                                    kaon_XbVsPhT_rec->Fill(xbj_DA_kaon_Rec, PhT_kaon_Rec);
                                    kaon_XbVsZ_rec->Fill(xbj_DA_kaon_Rec, z_DA_kaon_Rec);
                                    kaon_XbVsPhi_rec->Fill(reckaonPhi, xbj_DA_kaon_Rec);
                                    kaon_XbVsTheta_rec->Fill(reckaonTheta, xbj_DA_kaon_Rec);
                                    if(recpdg==pdg){
                                      RealkaonQ2->Fill(Q2_DA_kaon_Rec);
                                      RealkaonX->Fill(xbj_DA_kaon_Rec);
                                      RealkaonZ->Fill(z_DA_kaon_Rec);
                                      RealkaonPhT->Fill(PhT_kaon_Rec);
                                      RealkaonEta->Fill(reckaonEta);
                                      RealkaonPhi->Fill(reckaonPhi);
                                      RealkaonMom->Fill(mom_kaon_Rec);
                                    }
                                  }
                                  PIDgen2rec->Fill(kkaon, checkRecoID(recpdg));
                                }
                              }
                            }                            
                          }
                      }
                    scatElq_kaon.push_back(currentQ2kaon);
                    double y_DA_kaon = ((TMath::Tan(kaonPhiRad*0.5)) / (TMath::Tan(kaonPhiRad*0.5) + TMath::Tan(currentPhi*0.5)));
                    if(y_DA_kaon <= 0.95){
                        if(parentsIndex[i]<=5){
                            double Q2_DA_kaon = ((4*18*18*(1-y_DA_kaon)) / (TMath::Tan(currentPhi*0.5)*TMath::Tan(currentPhi*0.5)));
                            double xbj_DA_kaon = (Q2_DA_kaon) / (4*18*275*y_DA_kaon);
                            double z_DA_kaon = (ProtonBeam * kaon) / (ProtonBeam * photon_kaon);
                            TVector3 PhT_kaon_vec = trueKaonMom - ((trueKaonMom * currentQ2kaon) / currentQ2kaon.Mag())*currentQ2kaon.Unit();
                            double PhT_kaon = PhT_kaon_vec.Mag();
                            if(z_DA_kaon <= 1){
                              truekaonQ2->Fill(Q2_DA_kaon);
                              truekaon_xbj->Fill(xbj_DA_kaon);
                              truekaon_z->Fill(z_DA_kaon);
                              truekaon_PhT->Fill(PhT_kaon);
                              truekaon_mom->Fill(mom_kaon);
                              truekaonEta->Fill(kaonEta);
                              truekaonPhi->Fill(kaonPhi);
                              kaon_MomVsEta->Fill(kaonEta, mom_kaon);
                              kaon_PhTvsEta->Fill(kaonEta, PhT_kaon);
                              kaon_MomVsPhi->Fill(kaonPhi, mom_kaon);
                              kaon_PhTvsZ->Fill(z_DA_kaon, PhT_kaon);
                              kaon_MomVsTheta->Fill(kaonTheta, mom_kaon);
                              kaon_ThetaVsPhi->Fill(kaonTheta,kaonPhi);
                              kaon_PhTvsPhi->Fill(kaonPhi, PhT_kaon);
                              kaon_PhTvsTheta->Fill(kaonTheta, PhT_kaon);
                              //kaon_MomVsPhiVsTheta->Fill(kaonTheta, mom_kaon, kaonPhi);
                              kaon_PhTvsMom->Fill(PhT_kaon,mom_kaon);
                              kaon_ZvsMom->Fill(z_DA_kaon, mom_kaon);
                              kaon_ZvsPhi->Fill(kaonPhi,z_DA_kaon);
                              kaon_ZvsTheta->Fill(kaonTheta, z_DA_kaon);
                              kaon_XbVsMom->Fill(xbj_DA_kaon, mom_kaon);
                              kaon_XbVsPhT->Fill(xbj_DA_kaon, PhT_kaon);
                              kaon_XbVsZ->Fill(xbj_DA_kaon, z_DA_kaon);
                              kaon_XbVsPhi->Fill(kaonPhi, xbj_DA_kaon);
                              kaon_XbVsTheta->Fill(kaonTheta, xbj_DA_kaon);
                              //kaon_PhTvsZvsXb->Fill(xbj_DA_kaon, z_DA_kaon, PhT_kaon); 
                            }
                        }
                    
                    }
                  }   
                }
              // ELECTRON PDG   
              if(pdg == 11)
                {
                  for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
                    {
                      int recpdg = recPdg[j];
                      if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
                        {
                          if(parentsIndex[i]<=5){
                            PIDgen2rec->Fill(kelectron, checkRecoID(recpdg));
                          }
                        }
                    }
                }
              // PROTON
              if(pdg == 2212)
                {
                  TVector3 trueProtonMom(partMomX[i],partMomY[i],partMomZ[i]);
                  float mom_proton = trueProtonMom.Mag();
                  float E_proton = sqrt(mom_proton*mom_proton + 0.937*0.937);
                  TLorentzVector proton(E_proton, partMomX[i],partMomY[i],partMomZ[i]);
                  float protonEta = trueProtonMom.PseudoRapidity();

                  if(protonEta <= 3.5 && protonEta >= -3.5){
                    double protonPhiRad = trueProtonMom.Theta();
                    float protonPhi = protonPhiRad * (180.0 / TMath::Pi());
        
                    TLorentzVector scatElproton(currentMom, currentQ2proton.X(), currentQ2proton.Y(), currentQ2proton.Z());
                    TLorentzVector photon_proton = ElectronBeam - scatElproton;

                    for(unsigned int j=0; j<simuAssoc.GetSize(); j++)
                      {
                        int recpdg = std::abs(recPdg[j]);
                        if(simuAssoc[j] == i) // Find association index matching the index of the thrown particle we are looking at
                          {
                            TVector3 recProtMom(trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]); 
                            float mom_proton_Rec = recProtMom.Mag();
                            float E_proton_rec = sqrt(mom_proton_Rec*mom_proton_Rec + 0.937*0.937);
                            TLorentzVector proton_Rec(E_proton_rec, trackMomX[recoAssoc[j]],trackMomY[recoAssoc[j]],trackMomZ[recoAssoc[j]]);
                            float recprotEta = recProtMom.PseudoRapidity();

                            if(recprotEta <= 3.5 && recprotEta >= -3.5){
                              double recprPhi = recProtMom.Theta();
                              float recprotPhi = recprPhi * (180.0 / TMath::Pi());
                              TLorentzVector scatElproton_Rec(currentMom, currentQ2proton.X(), currentQ2proton.Y(), currentQ2proton.Z());
                              TLorentzVector photon_proton_Rec = ElectronBeam - scatElproton_Rec;

                              double y_DA_proton_Rec = ((TMath::Tan(recprPhi*0.5)) / (TMath::Tan(recprPhi*0.5) + TMath::Tan(currentPhi*0.5)));
                              if(y_DA_proton_Rec <= 0.95){
                                if(parentsIndex[i]<=5){
                                  double Q2_DA_proton_Rec = ((4*18*18*(1-y_DA_proton_Rec)) / (TMath::Tan(currentPhi*0.5)*TMath::Tan(currentPhi*0.5)));
                                  double xbj_DA_proton_Rec = (Q2_DA_proton_Rec) / (4*18*275*y_DA_proton_Rec);
                                  double z_DA_proton_Rec = (ProtonBeam * proton_Rec) / (ProtonBeam * photon_proton_Rec);
                                  TVector3 PhT_proton_vec_Rec = recProtMom - ((recProtMom * currentQ2proton) / currentQ2proton.Mag())*currentQ2proton.Unit();
                                  double PhT_proton_Rec = PhT_proton_vec_Rec.Mag();
                                  if(z_DA_proton_Rec <= 1){
                                    RecprotonQ2->Fill(Q2_DA_proton_Rec);
                                    Recproton_xbj->Fill(xbj_DA_proton_Rec);
                                    Recproton_z->Fill(z_DA_proton_Rec);
                                    Recproton_PhT->Fill(PhT_proton_Rec);
                                    Recproton_mom->Fill(mom_proton_Rec);
                                    RecprotonEta->Fill(recprotEta);
                                    RecprotonPhi->Fill(recprotPhi);
                                    proton_MomVsEta_rec->Fill(recprotEta, mom_proton_Rec);
                                    proton_MomVsPhi_rec->Fill(recprotPhi, mom_proton_Rec);
                                    if(recpdg==pdg){
                                      RealprotonQ2->Fill(Q2_DA_proton_Rec);
                                      RealprotonX->Fill(xbj_DA_proton_Rec);
                                      RealprotonZ->Fill(z_DA_proton_Rec);
                                      RealprotonPhT->Fill(PhT_proton_Rec);
                                      RealprotonEta->Fill(recprotEta);
                                      RealprotonPhi->Fill(recprotPhi);
                                      RealprotonMom->Fill(mom_proton_Rec);
                                    }
                                  }
                                  PIDgen2rec->Fill(kproton, checkRecoID(recpdg));
                                }
                              }
                            }
                          }
                      }

                    scatElq_proton.push_back(currentQ2proton);
                    double y_DA_proton = ((TMath::Tan(protonPhiRad*0.5)) / (TMath::Tan(protonPhiRad*0.5) + TMath::Tan(currentPhi*0.5)));
                    if(y_DA_proton <= 0.95){ 
                        if(parentsIndex[i]<=5){
                            double Q2_DA_proton = ((4*18*18*(1-y_DA_proton)) / (TMath::Tan(currentPhi*0.5)*TMath::Tan(currentPhi*0.5)));
                            double xbj_DA_proton = (Q2_DA_proton) / (4*18*275*y_DA_proton);
                            double z_DA_proton = (ProtonBeam * proton) / (ProtonBeam * photon_proton);
                            TVector3 PhT_proton_vec = trueProtonMom - ((trueProtonMom * currentQ2proton) / currentQ2proton.Mag())*currentQ2proton.Unit();
                            double PhT_proton = PhT_proton_vec.Mag();
                            if(z_DA_proton <= 1){
                            trueprotonQ2->Fill(Q2_DA_proton);
                            trueproton_xbj->Fill(xbj_DA_proton);
                            trueproton_z->Fill(z_DA_proton);
                            trueproton_PhT->Fill(PhT_proton);
                            trueproton_mom->Fill(mom_proton);
                            trueprotonEta->Fill(protonEta);
                            trueprotonPhi->Fill(protonPhi);
                            proton_MomVsEta->Fill(protonEta, mom_proton);
                            proton_MomVsPhi->Fill(protonPhi, mom_proton);
                                
                            }   
                        }
                    }
                  }      
                }
            }
        }
    } 
    /*
    for(int iPi = 1; iPi <= pion_count; iPi++){
      double eStat_uPDF = CalcolaErroreStatistico(val_uPDF);
      double eStat_dPDF = CalcolaErroreStatistico(val_dPDF);
      double eStat_glPDF = CalcolaErroreStatistico(val_glPDF);
      double estat_x = CalcolaErroreStatistico(val_x);
      uPDF->SetPointError(iPi, estat_x, eStat_uPDF);
      dPDF->SetPointError(iPi, estat_x, eStat_dPDF);
      glPDF->SetPointError(iPi, estat_x, 0.1*eStat_glPDF);
    }
    */

    // CANVAS DEL Q2 ______________________________________________________________________________________

    TCanvas *trueQ2 = new TCanvas("trueQ2", "Production with Q2", 800, 600);
    trueQ2->SetLogx();

    TH1D *truepionQ22 = (TH1D*)truepionQ2->Clone("truepionQ22");
    TH1D *truekaonQ22 = (TH1D*)truekaonQ2->Clone("truekaonQ22");
    TH1D *trueprotonQ22 = (TH1D*)trueprotonQ2->Clone("trueprotonQ22");

    
    double integral_pion11 = truepionQ2->Integral();
    double integral_kaon11 = truekaonQ2->Integral();
    double integral_proton11 = trueprotonQ2->Integral();
    double total11 = integral_kaon11 + integral_pion11 + integral_proton11;
    truepionQ2->Scale(1.0 / ( 30*total11));
    truekaonQ2->Scale(1.0 / ( 30*total11));
    trueprotonQ2->Scale(1.0 / ( 30*total11));
    /*
    int nBins = truepionQ2->GetNbinsX();
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = truepionQ2->GetBinContent(i) + truekaonQ2->GetBinContent(i) + trueprotonQ2->GetBinContent(i);

        if (sumBinContents > 0) {
            truepionQ2->SetBinContent(i, truepionQ2->GetBinContent(i) / sumBinContents);
            truekaonQ2->SetBinContent(i, truekaonQ2->GetBinContent(i) / sumBinContents);
            trueprotonQ2->SetBinContent(i, trueprotonQ2->GetBinContent(i) / sumBinContents);
        }
    }*/
    truepionQ2->SetStats(kFALSE);
    truepionQ2->SetLineColor(kRed);
    truepionQ2->SetTitle("MC Production of charged particles | 18x275 GeV");
    truepionQ2->GetXaxis()->SetTitle("Q^2 [GeV^2]");
    truepionQ2->GetYaxis()->SetTitle("Total fraction");
    truepionQ2->Draw("HIST");
    //truepionQ2->GetYaxis()->SetRangeUser(0, 1000);
    truekaonQ2->SetLineColor(kBlue);
    truekaonQ2->Draw("HIST SAME");
    trueprotonQ2->SetLineColor(kOrange);
    trueprotonQ2->Draw("HIST SAME");

    TLegend *legend2 = new TLegend(0.75, 0.78, 0.88, 0.88);
    legend2->AddEntry(truepionQ2, "Pions", "l");
    legend2->AddEntry(truekaonQ2, "Kaons", "l");
    legend2->AddEntry(trueprotonQ2, "Protons", "l");
    legend2->Draw();

    //trueQ2->SetTitle("Production of charged particles with Q2");
    trueQ2->Update();
    trueQ2->Write();

    TCanvas *RecQ2 = new TCanvas("RecQ2", "Reconstructed Q2", 800, 600);
    RecQ2->SetLogx();

    TH1D *RecpionQ22 = (TH1D*)RecpionQ2->Clone("RecpionQ22");
    TH1D *ReckaonQ22 = (TH1D*)ReckaonQ2->Clone("ReckaonQ22");
    TH1D *RecprotonQ22 = (TH1D*)RecprotonQ2->Clone("RecprotonQ22");
    
    double integral_pion12 = RecpionQ2->Integral();
    double integral_kaon12 = ReckaonQ2->Integral();
    double integral_proton12 = RecprotonQ2->Integral();
    double total12 = integral_kaon12 + integral_pion12 + integral_proton12;
    RecpionQ2->Scale(1.0 / ( 30*total12));
    ReckaonQ2->Scale(1.0 / ( 30*total12));
    RecprotonQ2->Scale(1.0 / ( 30*total12));
    /*
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = RecpionQ2->GetBinContent(i) + ReckaonQ2->GetBinContent(i) + RecprotonQ2->GetBinContent(i);

        if (sumBinContents > 0) {
            RecpionQ2->SetBinContent(i, RecpionQ2->GetBinContent(i) / sumBinContents);
            ReckaonQ2->SetBinContent(i, ReckaonQ2->GetBinContent(i) / sumBinContents);
            RecprotonQ2->SetBinContent(i, RecprotonQ2->GetBinContent(i) / sumBinContents);
        }
    }*/
    RecpionQ2->SetStats(kFALSE);
    RecpionQ2->SetLineColor(kRed);
    RecpionQ2->SetTitle("Reconstruction of charged particles | 18x275 GeV");
    RecpionQ2->GetXaxis()->SetTitle("Q^2 [GeV^2]");
    RecpionQ2->GetYaxis()->SetTitle("Total fraction");
    RecpionQ2->Draw("HIST");
    //RecpionQ2->GetYaxis()->SetRangeUser(0, 1000);
    ReckaonQ2->SetLineColor(kBlue);
    ReckaonQ2->Draw("HIST SAME");
    RecprotonQ2->SetLineColor(kOrange);
    RecprotonQ2->Draw("HIST SAME");

    legend2->Draw();

    RecQ2->Update();
    RecQ2->Write();

    TCanvas *S_Q2 = new TCanvas("Sensibility_Q2", "Sensibility of the Q2 measure", 800, 600);
    S_Q2->SetLogx();
    TH1D *Sinpione_Q2 = (TH1D*)RecpionQ22->Clone("Sinpione_Q2");
    Sinpione_Q2->Divide(truepionQ22);
    Sinpione_Q2->SetTitle("Efficiency reconstruction of charged particles | 18x275 GeV");
    Sinpione_Q2->GetXaxis()->SetTitle("Q^2");
    Sinpione_Q2->SetLineColor(kRed);
    Sinpione_Q2->GetYaxis()->SetTitle("Rec / MC events");
    Sinpione_Q2->SetStats(kFALSE);
    Sinpione_Q2->GetYaxis()->SetRangeUser(0, 2);
    Sinpione_Q2->Draw("HIST");
    TH1D *Sinkaon_Q2 = (TH1D*)ReckaonQ22->Clone("Sinkaon_Q2");
    Sinkaon_Q2->Divide(truekaonQ22);
    Sinkaon_Q2->SetLineColor(kBlue);
    Sinkaon_Q2->Draw("IHST SAME");
    TH1D *Sinproton_Q2 = (TH1D*)RecprotonQ22->Clone("Sinproton_Q2");
    Sinproton_Q2->Divide(trueprotonQ22);
    Sinproton_Q2->SetLineColor(kOrange);
    Sinproton_Q2->Draw("HIST SAME");

    legend2->Draw();

    S_Q2->Update();
    S_Q2->Write();


    // CANVAS PER x_Bj_____________________________________________________________________________________________-

    TCanvas *true_xbj = new TCanvas("true_xbj", "Productions over x_Bj", 800, 600);
    true_xbj->SetLogx();

    TH1D *truepion_xbj2 = (TH1D*)truepion_xbj->Clone("truepion_xbj2");
    TH1D *truekaon_xbj2 = (TH1D*)truekaon_xbj->Clone("truekaon_xbj2");
    TH1D *trueproton_xbj2 = (TH1D*)trueproton_xbj->Clone("trueproton_xbj2");

    double integral_pion21 = truepion_xbj->Integral();
    double integral_kaon21 = truekaon_xbj->Integral();
    double integral_proton21 = trueproton_xbj->Integral();
    double total21 = integral_kaon21 + integral_pion21 + integral_proton21;
    truepion_xbj->Scale(1.0 / ( 30*total21));
    truekaon_xbj->Scale(1.0 / ( 30*total21));
    trueproton_xbj->Scale(1.0 / ( 30*total21));
    /*
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = truepion_xbj->GetBinContent(i) + truekaon_xbj->GetBinContent(i) + trueproton_xbj->GetBinContent(i);

        if (sumBinContents > 0) {
            truepion_xbj->SetBinContent(i, truepion_xbj->GetBinContent(i) / sumBinContents);
            truekaon_xbj->SetBinContent(i, truekaon_xbj->GetBinContent(i) / sumBinContents);
            trueproton_xbj->SetBinContent(i, trueproton_xbj->GetBinContent(i) / sumBinContents);
        }
    }*/
    truepion_xbj->SetStats(kFALSE);
    truepion_xbj->SetLineColor(kRed);
    truepion_xbj->SetTitle("MC Production of charged particles | 18x275 GeV");
    truepion_xbj->GetXaxis()->SetTitle("x_B");
    truepion_xbj->GetYaxis()->SetTitle("Total fraction");
    truepion_xbj->Draw("HIST");
    //truepion_xbj->GetYaxis()->SetRangeUser(0, 1500);
    truekaon_xbj->SetLineColor(kBlue);
    truekaon_xbj->Draw("HIST SAME");
    trueproton_xbj->SetLineColor(kOrange);
    trueproton_xbj->Draw("HIST SAME");

    legend2->Draw();

    true_xbj->Update();
    true_xbj->Write();

    TCanvas *Rec_xbj = new TCanvas("Rec_xbj", "Productions over x_Bj", 800, 600);
    Rec_xbj->SetLogx();
    TH1D *Recpion_xbj2 = (TH1D*)Recpion_xbj->Clone("Recpion_xbj2");
    TH1D *Reckaon_xbj2 = (TH1D*)Reckaon_xbj->Clone("Reckaon_xbj2");
    TH1D *Recproton_xbj2 = (TH1D*)Recproton_xbj->Clone("Recproton_xbj2");
    
    double integral_pion22 = Recpion_xbj->Integral();
    double integral_kaon22 = Reckaon_xbj->Integral();
    double integral_proton22 = Recproton_xbj->Integral();
    double total22 = integral_kaon22 + integral_pion22 + integral_proton22;
    Recpion_xbj->Scale(1.0 / ( 30*total22));
    Reckaon_xbj->Scale(1.0 / ( 30*total22));
    Recproton_xbj->Scale(1.0 / ( 30*total22));
    /*
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = Recpion_xbj->GetBinContent(i) + Reckaon_xbj->GetBinContent(i) + Recproton_xbj->GetBinContent(i);

        if (sumBinContents > 0) {
            Recpion_xbj->SetBinContent(i, Recpion_xbj->GetBinContent(i) / sumBinContents);
            Reckaon_xbj->SetBinContent(i, Reckaon_xbj->GetBinContent(i) / sumBinContents);
            Recproton_xbj->SetBinContent(i, Recproton_xbj->GetBinContent(i) / sumBinContents);
        }
    }*/
    Recpion_xbj->SetStats(kFALSE);
    Recpion_xbj->SetLineColor(kRed);
    Recpion_xbj->SetTitle("Reconstruction of charged particles| 18x275 GeV");
    Recpion_xbj->GetXaxis()->SetTitle("x_B");
    Recpion_xbj->GetYaxis()->SetTitle("Total fraction");
    Recpion_xbj->Draw("HIST");
    //Recpion_xbj->GetYaxis()->SetRangeUser(0, 1500);
    Reckaon_xbj->SetLineColor(kBlue);
    Reckaon_xbj->Draw("HIST SAME");
    Recproton_xbj->SetLineColor(kOrange);
    Recproton_xbj->Draw("HIST SAME");

    legend2->Draw();

    Rec_xbj->Update();
    Rec_xbj->Write();

    TCanvas *S_xbj = new TCanvas("Sensibility_xbj", "Sensibility of the x_Bj measure", 800, 600);
    S_xbj->SetLogx();
    TH1D *Sinpione_xbj = (TH1D*)Recpion_xbj2->Clone("Sinpione_xbj");
    Sinpione_xbj->Divide(truepion_xbj2);
    Sinpione_xbj->SetStats(kFALSE);
    Sinpione_xbj->SetTitle("Efficiency reconstruction of charged particles | 18x275 GeV");
    Sinpione_xbj->GetYaxis()->SetRangeUser(0, 2);
    Sinpione_xbj->GetXaxis()->SetTitle("x_B");
    Sinpione_xbj->GetYaxis()->SetTitle("Rec / MC events");
    Sinpione_xbj->SetLineColor(kRed);
    Sinpione_xbj->Draw("HIST");
    TH1D *Sinkaon_xbj = (TH1D*)Reckaon_xbj2->Clone("Sinkaon_xbj");
    Sinkaon_xbj->Divide(truekaon_xbj2);
    Sinkaon_xbj->SetLineColor(kBlue);
    Sinkaon_xbj->Draw("HIST SAME");
    TH1D *Sinproton_xbj = (TH1D*)Recproton_xbj2->Clone("Sinproton_xbj");
    Sinproton_xbj->Divide(trueproton_xbj2);
    Sinproton_xbj->SetLineColor(kOrange);
    Sinproton_xbj->Draw("HIST SAME");

    legend2->Draw();

    S_xbj->Update();
    S_xbj->Write();

    // CANVAS PER z __________________________________________________________________________________________

    TCanvas *true_z = new TCanvas("true_z", "Productions over z", 800, 600);
    true_z->SetLogx();

    TH1D *truepion_z2 = (TH1D*)truepion_z->Clone("truepion_z2");
    TH1D *truekaon_z2 = (TH1D*)truekaon_z->Clone("truekaon_z2");
    TH1D *trueproton_z2 = (TH1D*)trueproton_z->Clone("trueproton_z2");
    
    double integral_pion31 = truepion_z->Integral();
    double integral_kaon31 = truekaon_z->Integral();
    double integral_proton31 = trueproton_z->Integral();
    double total31 = integral_kaon31 + integral_pion31 + integral_proton31;
    truepion_z->Scale(1.0 / ( 30*total31));
    truekaon_z->Scale(1.0 / ( 30*total31));
    trueproton_z->Scale(1.0 / ( 30*total31));
    /*
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = truepion_z->GetBinContent(i) + truekaon_z->GetBinContent(i) + trueproton_z->GetBinContent(i);

        if (sumBinContents > 0) {
            truepion_z->SetBinContent(i, truepion_z->GetBinContent(i) / sumBinContents);
            truekaon_z->SetBinContent(i, truekaon_z->GetBinContent(i) / sumBinContents);
            trueproton_z->SetBinContent(i, trueproton_z->GetBinContent(i) / sumBinContents);
        }
    }
    */
    truepion_z->SetStats(kFALSE);
    truepion_z->SetLineColor(kRed);
    truepion_z->SetTitle("MC Production of charged particles | 18x275 GeV");
    truepion_z->GetXaxis()->SetTitle("z");
    truepion_z->GetYaxis()->SetTitle("Total fraction");
    truepion_z->Draw("HIST");
    //Recpion_z->GetYaxis()->SetRangeUser(0, 1500);
    truekaon_z->SetLineColor(kBlue);
    truekaon_z->Draw("HIST SAME");
    trueproton_z->SetLineColor(kOrange);
    trueproton_z->Draw("HIST SAME");
  
    TLegend *legend = new TLegend(0.17, 0.7, 0.3, 0.88);
    legend->AddEntry(truepion_z, "Pions", "l");
    legend->AddEntry(truekaon_z, "Kaons", "l");
    legend->AddEntry(trueproton_z, "Protons", "l");
    legend2->Draw();

    true_z->Update();
    true_z->Write();

    TCanvas *Rec_z = new TCanvas("Rec_z", "Productions over z", 800, 600);
    Rec_z->SetLogx();
    TH1D *Recpion_z2 = (TH1D*)Recpion_z->Clone("Recpion_z2");
    TH1D *Reckaon_z2 = (TH1D*)Reckaon_z->Clone("Reckaon_z2");
    TH1D *Recproton_z2 = (TH1D*)Recproton_z->Clone("Recproton_z2");
    
    double integral_pion32 = Recpion_z->Integral();
    double integral_kaon32 = Reckaon_z->Integral();
    double integral_proton32 = Recproton_z->Integral();
    double total32 = integral_kaon32 + integral_pion32 + integral_proton32;
    Recpion_z->Scale(1.0 / ( 30*total32));
    Reckaon_z->Scale(1.0 / ( 30*total32));
    Recproton_z->Scale(1.0 / ( 30*total32));
    /*
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = Recpion_z->GetBinContent(i) + Reckaon_z->GetBinContent(i) + Recproton_z->GetBinContent(i);

        if (sumBinContents > 0) {
            Recpion_z->SetBinContent(i, Recpion_z->GetBinContent(i) / sumBinContents);
            Reckaon_z->SetBinContent(i, Reckaon_z->GetBinContent(i) / sumBinContents);
            Recproton_z->SetBinContent(i, Recproton_z->GetBinContent(i) / sumBinContents);
        }
    }  
    */
    Recpion_z->SetStats(kFALSE);
    Recpion_z->SetLineColor(kRed);
    Recpion_z->SetTitle("Reconstruction of charged particles | 18x275 GeV");
    Recpion_z->GetXaxis()->SetTitle("z");
    Recpion_z->GetYaxis()->SetTitle("Total fraction");
    Recpion_z->Draw("HIST");
    //Recpion_z->GetYaxis()->SetRangeUser(0, 1500);
    Reckaon_z->SetLineColor(kBlue);
    Reckaon_z->Draw("HIST SAME");
    Recproton_z->SetLineColor(kOrange);
    Recproton_z->Draw("HIST SAME");

    legend2->Draw();

    Rec_z->Update();
    Rec_z->Write();

    TCanvas *S_z = new TCanvas("Sensibility_z", "Sensibility of the z measure", 800, 600);
    S_z->SetLogx();
    TH1D *Sinpione_z = (TH1D*)Recpion_z2->Clone("Sinpione_z");
    Sinpione_z->Divide(truepion_z2);
    Sinpione_z->SetStats(kFALSE);
    Sinpione_z->SetTitle("Efficiency reconstruction of charged particles | 18x275 GeV");
    Sinpione_z->GetXaxis()->SetTitle("z");
    Sinpione_z->GetYaxis()->SetTitle("Rec / MC events");
    Sinpione_z->SetLineColor(kRed);
    Sinpione_z->GetYaxis()->SetRangeUser(0, 2);
    Sinpione_z->Draw("HIST");
    TH1D *Sinkaon_z = (TH1D*)Reckaon_z2->Clone("Sinkaon_z");
    Sinkaon_z->Divide(truekaon_z2);
    Sinkaon_z->SetLineColor(kBlue);
    Sinkaon_z->Draw("HIST SAME");
    TH1D *Sinproton_z = (TH1D*)Recproton_z2->Clone("Sinproton_z");
    Sinproton_z->Divide(trueproton_z2);
    Sinproton_z->SetLineColor(kOrange);
    Sinproton_z->Draw("HIST SAME");

    legend2->Draw();

    S_z->Update();
    S_z->Write();

    // CANVAS PER P_hT _______________________________________________________________________________________________________

    TCanvas *true_PhT = new TCanvas("true_PhT", "Production over P_hT;", 800, 600);
    true_PhT->SetLogx();

    TH1D *truepion_PhT2 = (TH1D*)truepion_PhT->Clone("truepion_PhT2");
    TH1D *truekaon_PhT2 = (TH1D*)truekaon_PhT->Clone("truekaon_PhT2");
    TH1D *trueproton_PhT2 = (TH1D*)trueproton_PhT->Clone("trueproton_PhT2");
    
    double integral_pion41 = truepion_PhT->Integral();
    double integral_kaon41 = truekaon_PhT->Integral();
    double integral_proton41 = trueproton_PhT->Integral();
    double total41 = integral_kaon41 + integral_pion41 + integral_proton41;
    truepion_PhT->Scale(1.0 / ( 30*total41));
    truekaon_PhT->Scale(1.0 / ( 30*total41));
    trueproton_PhT->Scale(1.0 / ( 30*total41));
    /*
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = truepion_PhT->GetBinContent(i) + truekaon_PhT->GetBinContent(i) + trueproton_PhT->GetBinContent(i);

        if (sumBinContents > 0) {
            truepion_PhT->SetBinContent(i, truepion_PhT->GetBinContent(i) / sumBinContents);
            truekaon_PhT->SetBinContent(i, truekaon_PhT->GetBinContent(i) / sumBinContents);
            trueproton_PhT->SetBinContent(i, trueproton_PhT->GetBinContent(i) / sumBinContents);
        }
    }
    */
    truepion_PhT->SetStats(kFALSE);
    truepion_PhT->SetLineColor(kRed);
    truepion_PhT->SetTitle("MC Production of charged particles | 18x275 GeV");
    truepion_PhT->GetXaxis()->SetTitle("P_hT [GeV]");
    truepion_PhT->GetYaxis()->SetTitle("Total fraction");
    //truepion_PhT->GetYaxis()->SetRangeUser(0, 0.01);
    truepion_PhT->Draw("HIST");
    truekaon_PhT->SetLineColor(kBlue);
    truekaon_PhT->Draw("HIST SAME");
    trueproton_PhT->SetLineColor(kOrange);
    trueproton_PhT->Draw("HIST SAME");

    legend2->Draw();

    true_PhT->Update();
    true_PhT->Write();

    TCanvas *Rec_PhT = new TCanvas("Rec_PhT", "Production over P_hT;", 800, 600);
    Rec_PhT->SetLogx();

    TH1D *Recpion_PhT2 = (TH1D*)Recpion_PhT->Clone("Recpion_PhT2");
    TH1D *Reckaon_PhT2 = (TH1D*)Reckaon_PhT->Clone("Reckaon_PhT2");
    TH1D *Recproton_PhT2 = (TH1D*)Recproton_PhT->Clone("Recproton_PhT2");

    double integral_pion42 = Recpion_PhT->Integral();
    double integral_kaon42 = Reckaon_PhT->Integral();
    double integral_proton42 = Recproton_PhT->Integral();
    double total42 = integral_kaon42 + integral_pion42 + integral_proton42;
    Recpion_PhT->Scale(1.0 / ( 30*total42));
    Reckaon_PhT->Scale(1.0 / ( 30*total42));
    Recproton_PhT->Scale(1.0 / ( 30*total42));
    /*
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = Recpion_PhT->GetBinContent(i) + Reckaon_PhT->GetBinContent(i) + Recproton_PhT->GetBinContent(i);

        if (sumBinContents > 0) {
            Recpion_PhT->SetBinContent(i, Recpion_PhT->GetBinContent(i) / sumBinContents);
            Reckaon_PhT->SetBinContent(i, Reckaon_PhT->GetBinContent(i) / sumBinContents);
            Recproton_PhT->SetBinContent(i, Recproton_PhT->GetBinContent(i) / sumBinContents);
        }
    }
    */
    Recpion_PhT->SetStats(kFALSE);
    Recpion_PhT->SetLineColor(kRed);
    Recpion_PhT->SetTitle("Recontruction of charged particles | 18x275 GeV");
    Recpion_PhT->GetXaxis()->SetTitle("P_hT [GeV]");
    Recpion_PhT->GetYaxis()->SetTitle("Total fraction"); 
    //Recpion_PhT->GetYaxis()->SetRangeUser(0, 0.01);
    Recpion_PhT->Draw("HIST");
    Reckaon_PhT->SetLineColor(kBlue);
    Reckaon_PhT->Draw("HIST SAME");
    Recproton_PhT->SetLineColor(kOrange);
    Recproton_PhT->Draw("HIST SAME");

    legend2->Draw();

    Rec_PhT->Update();
    Rec_PhT->Write();

    TCanvas *S_PhT = new TCanvas("Sensibility_PhT", "Sensibility of the P_hT measure", 800, 600);
    S_PhT->SetLogx();
    TH1D *Sinpione_PhT = (TH1D*)Recpion_PhT2->Clone("Sinpione_PhT");
    Sinpione_PhT->Divide(truepion_PhT2);
    Sinpione_PhT->SetStats(kFALSE);
    Sinpione_PhT->SetLineColor(kRed);
    Sinpione_PhT->SetTitle("Efficiency reconstruction of charged particles | 18x275 GeV");
    Sinpione_PhT->GetXaxis()->SetTitle("P_hT [GeV]");
    Sinpione_PhT->GetYaxis()->SetTitle("Rec / MC events");
    Sinpione_PhT->GetYaxis()->SetRangeUser(0, 2);
    Sinpione_PhT->Draw("HIST");
    TH1D *Sinkaon_PhT = (TH1D*)Reckaon_PhT2->Clone("Sinkaon_PhT");
    Sinkaon_PhT->Divide(truekaon_PhT2);
    Sinkaon_PhT->SetLineColor(kBlue);
    Sinkaon_PhT->Draw("HIST SAME");
    TH1D *Sinproton_PhT = (TH1D*)Recproton_PhT2->Clone("Sinproton_PhT");
    Sinproton_PhT->Divide(trueproton_PhT2);
    Sinproton_PhT->SetLineColor(kOrange);
    Sinproton_PhT->Draw("HIST SAME");

    legend2->Draw();

    S_PhT->Update();
    S_PhT->Write();

    //  CANVAS PER ETA ________________________________________________________________________________________________________

    TCanvas *trueEtaDist = new TCanvas("trueEtaDist", "Eta Distributions", 800, 600);
    TH1D *truepionEta2 = (TH1D*)truepionEta->Clone("truepionEta2");
    TH1D *truekaonEta2 = (TH1D*)truekaonEta->Clone("truekaonEta2");
    TH1D *trueprotonEta2 = (TH1D*)trueprotonEta->Clone("trueprotonEta2");

    double integral_pion51 = truepionEta->Integral();
    double integral_kaon51 = truekaonEta->Integral();
    double integral_proton51 = trueprotonEta->Integral();
    double total51 = integral_kaon51 + integral_pion51 + integral_proton51;
    truepionEta->Scale(1.0 / ( 30*total51));
    truekaonEta->Scale(1.0 / ( 30*total51));
    trueprotonEta->Scale(1.0 / ( 30*total51));
     /*for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = truepionEta->GetBinContent(i) + truekaonEta->GetBinContent(i) + trueprotonEta->GetBinContent(i);

        if (sumBinContents > 0) {
            truepionEta->SetBinContent(i, truepionEta->GetBinContent(i) / sumBinContents);
            truekaonEta->SetBinContent(i, truekaonEta->GetBinContent(i) / sumBinContents);
            trueprotonEta->SetBinContent(i, trueprotonEta->GetBinContent(i) / sumBinContents);
        }
    }*/
    truepionEta->SetLineColor(kRed);
    truepionEta->SetStats(kFALSE);
    truepionEta->SetTitle("MC Production of charged particles | 18x275 GeV");
    truepionEta->GetXaxis()->SetTitle("Eta");
    truepionEta->GetYaxis()->SetTitle("Total fraction");
    //truepionEta->GetYaxis()->SetRangeUser(0, 0.08);

    truepionEta->Draw("HIST");
    truekaonEta->SetLineColor(kBlue);
    truekaonEta->Draw("HIST SAME");
    trueprotonEta->SetLineColor(kOrange);
    trueprotonEta->Draw("HIST SAME");

    legend2->Draw();
    
    trueEtaDist->Update();
    trueEtaDist->Write();

    TCanvas *RecEtaDist = new TCanvas("RecEtaDist", "Eta Distributions", 800, 600);

    TH1D *RecpionEta2 = (TH1D*)RecpionEta->Clone("RecpionEta2");
    TH1D *ReckaonEta2 = (TH1D*)ReckaonEta->Clone("ReckaonEta2");
    TH1D *RecprotonEta2 = (TH1D*)RecprotonEta->Clone("RecprotonEta2");

    double integral_pion52 = RecpionEta->Integral();
    double integral_kaon52 = ReckaonEta->Integral();
    double integral_proton52 = RecprotonEta->Integral();
    double total52 = integral_kaon52 + integral_pion52 + integral_proton52;
    RecpionEta->Scale(1.0 / ( 30*total52));
    ReckaonEta->Scale(1.0 / ( 30*total52));
    RecprotonEta->Scale(1.0 / ( 30*total52));
    /*
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = RecpionEta->GetBinContent(i) + ReckaonEta->GetBinContent(i) + RecprotonEta->GetBinContent(i);

        if (sumBinContents > 0) {
            RecpionEta->SetBinContent(i, RecpionEta->GetBinContent(i) / sumBinContents);
            ReckaonEta->SetBinContent(i, ReckaonEta->GetBinContent(i) / sumBinContents);
            RecprotonEta->SetBinContent(i, RecprotonEta->GetBinContent(i) / sumBinContents);
        }
    }*/
    RecpionEta->SetLineColor(kRed);
    RecpionEta->SetStats(kFALSE);
    RecpionEta->SetTitle("Reconstruction of charged particles | 18x275 GeV");
    RecpionEta->GetXaxis()->SetTitle("Eta");
    RecpionEta->GetYaxis()->SetTitle("Total fraction");
    RecpionEta->Draw("HIST");
    ReckaonEta->SetLineColor(kBlue);
    ReckaonEta->Draw("HIST SAME");
    RecprotonEta->SetLineColor(kOrange);
    RecprotonEta->Draw("HIST SAME");

    // Creazione della leggenda
    legend2->Draw();
    
    RecEtaDist->Update();
    RecEtaDist->Write();

    TCanvas *S_Eta = new TCanvas("Sensibility_Eta", "Sensibility of the Eta measure", 800, 600);
    TH1D *Sinpione_Eta = (TH1D*)RecpionEta2->Clone("Sinpione_Eta");
    Sinpione_Eta->Divide(truepionEta2);
    Sinpione_Eta->SetStats(kFALSE);
    Sinpione_Eta->SetTitle("Efficiency reconstruction of charged particles | 18x275 GeV");
    Sinpione_Eta->GetXaxis()->SetTitle("Eta");
    Sinpione_Eta->GetYaxis()->SetTitle("Rec / MC events");
    Sinpione_Eta->SetLineColor(kRed);
    Sinpione_Eta->GetYaxis()->SetRangeUser(0, 2);
    Sinpione_Eta->Draw("HIST");
    TH1D *Sinkaon_Eta = (TH1D*)ReckaonEta2->Clone("Sinkaon_Eta");
    Sinkaon_Eta->Divide(truekaonEta2);
    Sinkaon_Eta->SetLineColor(kBlue);
    Sinkaon_Eta->Draw("HIST SAME");
    TH1D *Sinproton_Eta = (TH1D*)RecprotonEta2->Clone("Sinproton_Eta");
    Sinproton_Eta->Divide(trueprotonEta2);
    Sinproton_Eta->SetLineColor(kOrange);
    Sinproton_Eta->Draw("HIST SAME");

    legend2->Draw();

    S_Eta->Update();
    S_Eta->Write();

    // CANVAS PER PHI ___________________________________________________________________________________________________________________-
        
    TCanvas *truePhiDist = new TCanvas("truePhiDist", "Phi Distributions of generated particles", 800, 600);
    TH1D *truepionPhi2 = (TH1D*)truepionPhi->Clone("truepionPhi2");
    TH1D *truekaonPhi2 = (TH1D*)truekaonPhi->Clone("truekaonPhi2");
    TH1D *trueprotonPhi2 = (TH1D*)trueprotonPhi->Clone("trueprotonPhi2");

    double integral_pion61 = truepionPhi->Integral();
    double integral_kaon61 = truekaonPhi->Integral();
    double integral_proton61 = trueprotonPhi->Integral();
    double total61 = integral_kaon61 + integral_pion61 + integral_proton61;
    truepionPhi->Scale(1.0 / ( 30*total61));
    truekaonPhi->Scale(1.0 / ( 30*total61));
    trueprotonPhi->Scale(1.0 / ( 30*total61));
    /*
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = truepionPhi->GetBinContent(i) + truekaonPhi->GetBinContent(i) + trueprotonPhi->GetBinContent(i);

        if (sumBinContents > 0) {
            truepionPhi->SetBinContent(i, truepionPhi->GetBinContent(i) / sumBinContents);
            truekaonPhi->SetBinContent(i, truekaonPhi->GetBinContent(i) / sumBinContents);
            trueprotonPhi->SetBinContent(i, trueprotonPhi->GetBinContent(i) / sumBinContents);
        }
    }*/
    truepionPhi->SetLineColor(kRed);
    truepionPhi->SetStats(kFALSE);
    truepionPhi->SetTitle("MC Production of charged particles | 18x275 GeV");
    truepionPhi->GetXaxis()->SetTitle("Phi");
    truepionPhi->GetYaxis()->SetTitle("Total fraction");
   //truepionPhi->GetYaxis()->SetRangeUser(0, 1);
    truepionPhi->Draw("HIST");
    truekaonPhi->SetLineColor(kBlue);
    truekaonPhi->Draw("HIST SAME");
    trueprotonPhi->SetLineColor(kOrange);
    trueprotonPhi->Draw("HIST SAME");

    // Creazione della leggenda
    legend2->Draw();

    truePhiDist->Update();
    truePhiDist->Write();
    
    TCanvas *RecPhiDist = new TCanvas("RecPhiDist", "Phi Distributions of reconstructed particles", 800, 600);
    TH1D *RecpionPhi2 = (TH1D*)RecpionPhi->Clone("RecpionPhi2");
    TH1D *ReckaonPhi2 = (TH1D*)ReckaonPhi->Clone("ReckaonPhi2");
    TH1D *RecprotonPhi2 = (TH1D*)RecprotonPhi->Clone("RecprotonPhi2");

    double integral_pion62 = RecpionPhi->Integral();
    double integral_kaon62 = ReckaonPhi->Integral();
    double integral_proton62 = RecprotonPhi->Integral();
    double total62 = integral_kaon62 + integral_pion62 + integral_proton62;
    RecpionPhi->Scale(1.0 / ( 30*total62));
    ReckaonPhi->Scale(1.0 / ( 30*total62));
    RecprotonPhi->Scale(1.0 / ( 30*total62));
    /*
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = RecpionPhi->GetBinContent(i) + ReckaonPhi->GetBinContent(i) + RecprotonPhi->GetBinContent(i);

        if (sumBinContents > 0) {
            RecpionPhi->SetBinContent(i, RecpionPhi->GetBinContent(i) / sumBinContents);
            ReckaonPhi->SetBinContent(i, ReckaonPhi->GetBinContent(i) / sumBinContents);
            RecprotonPhi->SetBinContent(i, RecprotonPhi->GetBinContent(i) / sumBinContents);
        }
    } */
    RecpionPhi->SetLineColor(kRed);
    RecpionPhi->SetStats(kFALSE);
    RecpionPhi->SetTitle("Reconstruction of charged particles | 18x275 GeV");
    RecpionPhi->GetXaxis()->SetTitle("Phi");
    RecpionPhi->GetYaxis()->SetTitle("Total fraction");
    RecpionPhi->Draw("HIST");
    ReckaonPhi->SetLineColor(kBlue);
    ReckaonPhi->Draw("HIST SAME");
    RecprotonPhi->SetLineColor(kOrange);
    RecprotonPhi->Draw("HIST SAME");

    // Creazione della leggenda

    legend2->Draw();

    RecPhiDist->Update();
    RecPhiDist->Write();
    
    TCanvas *S_Phi = new TCanvas("Sensibility_Phi", "Sensibility of the x_Bj measure", 800, 600);
    TH1D *Sinpione_Phi = (TH1D*)RecpionPhi2->Clone("Sinpione_Phi");
    Sinpione_Phi->Divide(truepionPhi2);
    Sinpione_Phi->SetStats(kFALSE);
    Sinpione_Phi->SetTitle("Efficiency reconstruction of charged particles | 18x275 GeV");
    Sinpione_Phi->GetXaxis()->SetTitle("Phi");
    Sinpione_Phi->GetYaxis()->SetTitle("Rec / MC events");
    Sinpione_Phi->SetLineColor(kRed);
    Sinpione_Phi->Draw("HIST");
    Sinpione_Phi->GetYaxis()->SetRangeUser(0, 2);
    TH1D *Sinkaon_Phi = (TH1D*)ReckaonPhi2->Clone("Sinkaon_Phi");
    Sinkaon_Phi->Divide(truekaonPhi2);
    Sinkaon_Phi->SetLineColor(kBlue);
    Sinkaon_Phi->Draw("HIST SAME");
    TH1D *Sinproton_Phi = (TH1D*)RecprotonPhi2->Clone("Sinproton_Phi");
    Sinproton_Phi->Divide(trueprotonPhi2);
    Sinproton_Phi->SetLineColor(kOrange);
    Sinproton_Phi->Draw("HIST SAME");

    legend2->Draw();

    S_Phi->Update();
    S_Phi->Write();

    // CANVAS PER IL MOMENTO ______________________________________________________________________________________________

    TCanvas *true_Mom = new TCanvas("True_Mom", "Production of charged particles with the Momentum", 800, 600);
    true_Mom->SetLogx();
    TH1D *truepion_mom2 = (TH1D*)truepion_mom->Clone("truepion_mom2");
    TH1D *truekaon_mom2 = (TH1D*)truekaon_mom->Clone("truekaon_mom2");
    TH1D *trueproton_mom2 = (TH1D*)trueproton_mom->Clone("trueproton_mom2");

    double integral_pion71 = truepion_mom->Integral();
    double integral_kaon71 = truekaon_mom->Integral();
    double integral_proton71 = trueproton_mom->Integral();
    double total71 = integral_kaon71 + integral_pion71 + integral_proton71;
    truepion_mom->Scale(1.0 / ( 30*total71));
    truekaon_mom->Scale(1.0 / ( 30*total71));
    trueproton_mom->Scale(1.0 / ( 30*total71));
   /*
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = truepion_mom->GetBinContent(i) + truekaon_mom->GetBinContent(i) + trueproton_mom->GetBinContent(i);

        if (sumBinContents > 0) {
            truepion_mom->SetBinContent(i, truepion_mom->GetBinContent(i) / sumBinContents);
            truekaon_mom->SetBinContent(i, truekaon_mom->GetBinContent(i) / sumBinContents);
            trueproton_mom->SetBinContent(i, trueproton_mom->GetBinContent(i) / sumBinContents);
        }
    }*/
    truepion_mom->SetLineColor(kRed);
    truepion_mom->SetTitle("MC production of charged particles | 18x275");
    truepion_mom->SetStats(kFALSE);
    truepion_mom->GetXaxis()->SetTitle("P_h [GeV]");
    truepion_mom->GetYaxis()->SetTitle("Total fraction");
    truepion_mom->Draw("HIST");
    truekaon_mom->SetLineColor(kBlue);
    truekaon_mom->Draw("HIST SAME");
    trueproton_mom->SetLineColor(kOrange);
    trueproton_mom->Draw("HIST SAME");

    legend2->Draw();

    true_Mom->Update();
    true_Mom->Write();

    TCanvas *Rec_Mom = new TCanvas("Rec_Mom", "Production of charged particles with the Momentum", 800, 600);
    Rec_Mom->SetLogx();
    TH1D *Recpion_mom2 = (TH1D*)Recpion_mom->Clone("Recpion_mom2");
    TH1D *Reckaon_mom2 = (TH1D*)Reckaon_mom->Clone("Reckaon_mom2");
    TH1D *Recproton_mom2 = (TH1D*)Recproton_mom->Clone("Recproton_mom2");

    double integral_pion72 = Recpion_mom->Integral();
    double integral_kaon72 = Reckaon_mom->Integral();
    double integral_proton72 = Recproton_mom->Integral();
    double total72 = integral_kaon72 + integral_pion72 + integral_proton72;
    Recpion_mom->Scale(1.0 / ( 30*total72));
    Reckaon_mom->Scale(1.0 / ( 30*total72));
    Recproton_mom->Scale(1.0 / ( 30*total72));
    /*
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = Recpion_mom->GetBinContent(i) + Reckaon_mom->GetBinContent(i) + Recproton_mom->GetBinContent(i);

        if (sumBinContents > 0) {
            Recpion_mom->SetBinContent(i, Recpion_mom->GetBinContent(i) / sumBinContents);
            Reckaon_mom->SetBinContent(i, Reckaon_mom->GetBinContent(i) / sumBinContents);
            Recproton_mom->SetBinContent(i, Recproton_mom->GetBinContent(i) / sumBinContents);
        }
    }
    */
    Recpion_mom->SetLineColor(kRed);
    Recpion_mom->SetTitle("Reconstruction of charged particles | 18x275");
    Recpion_mom->GetXaxis()->SetTitle("P_h [GeV]");
    Recpion_mom->SetStats(kFALSE);
    Recpion_mom->GetYaxis()->SetTitle("Total fraction");
    Recpion_mom->Draw("HIST");
    Reckaon_mom->SetLineColor(kBlue);
    Reckaon_mom->Draw("HIST SAME");
    Recproton_mom->SetLineColor(kOrange);
    Recproton_mom->Draw("HIST SAME");

    legend2->Draw();

    Rec_Mom->Update();
    Rec_Mom->Write();

    TCanvas *S_Mom = new TCanvas("Sensibility_mom", "Sensibility of the x_Bj measure", 800, 600);
    TH1D *Sinpione_mom = (TH1D*)Recpion_mom2->Clone("Sinpione_mom");
    S_Mom->SetLogx();
    Sinpione_mom->Divide(truepion_mom2);
    Sinpione_mom->SetStats(kFALSE);
    Sinpione_mom->SetTitle("Efficiency reconstruction of charged particles | 18x275 GeV");
    Sinpione_mom->SetLineColor(kRed);
    Sinpione_mom->GetYaxis()->SetTitle("Rec / MC events");
    Sinpione_mom->Draw("HIST");
    Sinpione_mom->GetYaxis()->SetRangeUser(0, 2);
    TH1D *Sinkaon_mom = (TH1D*)Reckaon_mom2->Clone("Sinkaon_mom");
    Sinkaon_mom->Divide(truekaon_mom2);
    Sinkaon_mom->SetLineColor(kBlue);
    Sinkaon_mom->Draw("HIST SAME");
    TH1D *Sinproton_mom = (TH1D*)Recproton_mom2->Clone("Sinproton_mom");
    Sinproton_mom->Divide(trueproton_mom2);
    Sinproton_mom->SetLineColor(kOrange);
    Sinproton_mom->Draw("HIST SAME");

    legend2->Draw();

    S_Mom->Update();
    S_Mom->Write();

    // PER IL GRAFICO IN 3D
    //gStyle->SetPalette(kRainBow);  // Usa la palette dei colori "RainBow"
    gStyle->SetOptStat(0);         // Disattiva le informazioni statistiche
    // Creazione di una nuova canvas
    TCanvas *c1 = new TCanvas("pion_XvsZvsPhT", "3D Histogram", 800, 600);
    c1->SetLogx();  
    c1->SetLogy();  
    c1->SetLogz();  
    // Disegna l'istogramma 3D con i colori sui bin
    pion_PhTvsZvsXb->Draw("BOX2Z");
    // Aggiorna il canvas per visualizzare tutto correttamente
    c1->Update();
    c1->Write();
    TCanvas *c2 = new TCanvas("pion_ThetavsMomvsPhi", "3D Histogram", 800, 600);
    //c1->SetLogx();  
    c2->SetLogy();  
    //c1->SetLogz();  
    // Disegna l'istogramma 3D con i colori sui bin
    pion_MomVsPhiVsTheta->Draw("BOX2Z");
    // Aggiorna il canvas per visualizzare tutto correttamente
    c2->Update();
    c2->Write();


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

    TCanvas *c4 = new TCanvas("c2", "Proiezione 2D", 800, 600);
    h2->Draw("COLZ"); 
    h3->Draw("COLZ SAME");
    c4->Update();
    c4->Write();

    /*
    TCanvas* c5 = new TCanvas("PDF1", "3D PDF Plot", 800, 600);
    c5->SetLogx();  // scala logaritmica per x
    c5->SetLogy();  // scala logaritmica per Q^2

    hPDF->Draw("COLZ");
    c5->Update();
    c5->Write();
    
    // PDF
    TCanvas* pdff = new TCanvas("PDF2", "PDF vs x", 800, 600);
    pdff->SetLogx();
    //pdff->SetLogy();
    uPDF->SetLineColor(kBlue); uPDF->SetTitle("u(x)"); uPDF->SetMarkerStyle(20); 
    uPDF->GetXaxis()->SetRangeUser(1e-2, 1); uPDF->GetYaxis()->SetRangeUser(0, 1);
    dPDF->SetLineColor(kRed); dPDF->SetTitle("d(x)");  dPDF->SetMarkerStyle(21); 
    dPDF->GetXaxis()->SetRangeUser(1e-2, 1); dPDF->GetYaxis()->SetRangeUser(0, 1);
    glPDF->SetLineColor(kGreen); glPDF->SetTitle("g(x)");  glPDF->SetMarkerStyle(22); 
    glPDF->GetXaxis()->SetRangeUser(1e-2, 1); glPDF->GetYaxis()->SetRangeUser(0, 1);

    uPDF->Draw("AP");  // Disegna il grafico della PDF per il quark up
    dPDF->Draw("P SAME"); // Disegna gli altri sullo stesso plot
    glPDF->Draw("P SAME");

    // Imposta i titoli degli assi e una leggenda
    uPDF->GetXaxis()->SetTitle("x");
    uPDF->GetYaxis()->SetTitle("x f(x, Q^2)");
    //pdff->BuildLegend();
    TLegend *legend8 = new TLegend(0.75, 0.78, 0.88, 0.88);
    legend8->AddEntry(uPDF, "u", "l");
    legend8->AddEntry(dPDF, "d", "l");
    legend8->AddEntry(glPDF, "g/10", "l");
    legend8->Draw();
    pdff->Update();
    pdff->Write();
    gPad->Modified(); 
    */

    std::cout << "particelle generate: " << count << std::endl;
    std::cout << "particelle generate nel range di rapidita': " << countEta << std::endl;
    std::cout << "the acceptance is: " << countEta / count << std::endl;
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

  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1933.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1932.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1931.eicrecon.tree.edm4eic.root";
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
  
  const char* inputFile1 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1927.eicrecon.tree.edm4eic.root";
  const char* inputFile2 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1926.eicrecon.tree.edm4eic.root";
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1925.eicrecon.tree.edm4eic.root";
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
  const char* inputFile3 = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.1887.eicrecon.tree.edm4eic.root";
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



