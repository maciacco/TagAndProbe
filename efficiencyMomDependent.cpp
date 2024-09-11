// based on https://github.com/RobertaCarla/tag-and-probe-studies

#include "Config.h"

using namespace std;

int dimension(TH1F *hist, int nbin_hist, float &min)
{
  int num;
  float width= hist->GetBinWidth(1);
  float bin_ledge[nbin_hist];
  for(int i=1; i<=nbin_hist; i++){
    bin_ledge[i]= hist->GetBinLowEdge(i);
    if( abs(bin_ledge[i]-0.)<width/2){
      min=i;
      num= min-1;
      break;
    }
  }
  return num;
}

TH1F* rebin(TH1F *g_hist, const char* name = "h_rebin")
{
  TH1F* h_rebin= new TH1F(name, "", kNbinsx, k_xbins);
  for(int i=1; i<=kNbinsx; i++){
    int xmin= g_hist->GetXaxis()->FindBin(k_xbins[i-1]+epsilon);
    int xmax= g_hist->GetXaxis()->FindBin(k_xbins[i]-epsilon);
    double integral= g_hist->Integral(xmin, xmax);
    h_rebin->SetBinContent(i, integral);
  }
  return h_rebin;
}

void efficiencyMomDependent(const char* inFileName = "AnalysisResults", const char* outFileName = "efficiencyMomDependentMC", const char* dir = "efficiency-q-a"){
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TFile fAnalysis(Form("%s.root", inFileName));
  TFile efficiencyMomDependent(Form("%s.root", outFileName), "RECREATE");

  TH3F *piRec = (TH3F*)fAnalysis.Get(Form("%s/piRec", dir));
  TH2F* piRec_xy =(TH2F*)piRec->Project3D("yx");
  TH1F *Decays = (TH1F*)piRec->ProjectionY("Decays",1,1,1,900);
  TH1F *ITS_TPC = (TH1F*)piRec->ProjectionY("ITS_TPC",6,6,1,900);

  int nbin = Decays->GetNbinsX();
  float min;
  int n = dimension(Decays, nbin, min);

  float maxX= Decays->GetXaxis()->GetBinUpEdge(nbin);
  TH1F *pi_pos_d = new TH1F("pi_pos_d", ";#it{p}_{T} (GeV/#it{c}); Entries", n, 0., maxX);
  TH1F *pi_neg_d = new TH1F("pi_neg_d", ";#it{p}_{T} (GeV/#it{c});Entries", n, 0., maxX);
  TH1F *pi_pos_IT = new TH1F("pi_pos_ITS_TPC", ";#it{p}_{T} (GeV/#it{c});Entries", n, 0., maxX);
  TH1F *pi_neg_IT = new TH1F("pi_neg_ITS_TPC", ";#it{p}_{T} (GeV/#it{c});Entries", n, 0., maxX);

  pi_pos_d->SetBinContent(1, 0);
  pi_neg_d->SetBinContent(1, 0);
  for(int i{2}; i <= n; i++){
    pi_pos_d->SetBinContent(i, Decays->GetBinContent(min - 1 + i));
    pi_neg_d->SetBinContent(i, Decays->GetBinContent(min - i));

    pi_pos_IT->SetBinContent(i, ITS_TPC->GetBinContent(min - 1 + i));
    pi_neg_IT->SetBinContent(i, ITS_TPC->GetBinContent(min - i));
  }

  TH1F* pi_neg_d_rebin = rebin(pi_neg_d, "pi_neg_d_rebin");
  TH1F* pi_neg_IT_rebin= rebin(pi_neg_IT, "pi_neg_IT_rebin");
  TH1F* pi_pos_d_rebin = rebin(pi_pos_d, "pi_pos_d_rebin");
  TH1F* pi_pos_IT_rebin= rebin(pi_pos_IT, "pi_pos_IT_rebin");

  TH1F *gE_glob_neg_rebin = new TH1F(*pi_neg_IT_rebin);
  gE_glob_neg_rebin->Divide(pi_neg_IT_rebin, pi_neg_d_rebin, 1., 1., "B");
  gE_glob_neg_rebin->SetName("EffPiMinusGlo");
  gE_glob_neg_rebin->SetTitle(";#it{p}_{T} (GeV/#it{c});#epsilon_{glo}(#pi^{-})");

  TH1F *gE_glob_pos_rebin = new TH1F(*pi_pos_IT_rebin);
  gE_glob_pos_rebin->Divide(pi_pos_IT_rebin, pi_pos_d_rebin, 1., 1., "B");
  gE_glob_pos_rebin->SetName("EffPiPlusGlo");
  gE_glob_pos_rebin->SetTitle(";#it{p}_{T} (GeV/#it{c});#epsilon_{glo}(#pi^{+})");

  gE_glob_neg_rebin->Write();
  gE_glob_pos_rebin->Write();
  efficiencyMomDependent.Close();
}