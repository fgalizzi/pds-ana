#include "../classe.hpp"
//*********************************************
void cla::Avg_Muon(){
//*********************************************
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); 

  std::vector<double> avg_mu;
  std::vector<std::vector<double>> sel_wf, mu_wf;

  read();
  
  TH1D* h_p = AllFpromptHisto(wfs, int_low, int_up, int_prompt);
  SelPDE_WF(wfs, sel_wf, prepulse_ticks, int_prompt, sat_low, amp_low, amp_up, bsl);  

  TH2D* hI = BuildChargeFpromptHisto(sel_wf, mu_wf, int_low, int_up,
      int_prompt, f_prompt);
  
  double ymin, ymax;
  min_max_element(mu_wf, ymin, ymax);

  TH2D* h2 = new TH2D("h2", Form("%s;%s;%s", "Persistence", "Ticks", "ADC Counts"),
      memorydepth/2, 0., memorydepth, 
      120, ymin, ymax);
  
  for (int i=0; i< mu_wf.size(); i++) 
    for (int j=0; j<memorydepth; j=j+2) h2->Fill(j, mu_wf[i][j]);

  TH1D* h_pr = hI->ProjectionY("h_pr");
  
  TCanvas *c1 = new TCanvas("c1","c1",0,0,500,450);  
  c1->cd();
  h2->Draw("COLZ");
  c1->SetLogz();
  c1->Modified(); c1->Update();


  TCanvas *c2 = new TCanvas("c2","c2",500,0,500,450);
  c2->cd();
  hI->Draw();
  c2->Modified(); c2->Update();

  avgWF(mu_wf, avg_mu);

  if(print==true) {
    VecDouble_in_Binary("Muon.dat", avg_mu);
    print = false;
  }

  TGraph *g1 = new TGraph(avg_mu.size(), &avg_mu[0]);
  //g1->GetXaxis()->SetTitle("Time [#mus]");
  g1->GetYaxis()->SetTitle("ADC counts");
  g1->SetLineColor(2);
  TCanvas *c3 = new TCanvas("c3","c3",0,550,500,450);
  c3->cd();
  g1->Draw();
  c3->Modified(); c3->Update();

  TCanvas *c4 = new TCanvas("c4","c4",500,550,500,450);
  c4->cd();
  h_p->SetLineColor(kRed);
  h_p->Draw();
  h_pr->Draw("SAME");
  c4->Modified(); c4->Update();
}
