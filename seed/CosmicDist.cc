#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
using namespace std;

#include "TApplication.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "TRandom.h"
#include "Setting.h"

const int NCanvas=2;

class CosmicDist
{
 public:
         CosmicDist();
        ~CosmicDist();
  void MakeHist(string ofname);
  void Loop();
  bool Draw(); 
  void SaveCanvas(string PdfFileName); 
  void SetFileName( string filename )  { ascii_filename = filename; }
  void OutputFileOpen();
  void WriteALine(TVector3 pos, TVector3 dir, int evID, int trID=1, int pareID=0, double weight=1.0);
  void SetDistType( int type )  { pos_type = type; }
  void SetMaxEvent( int N )  { ENumMax = N; }
  void SetDiameter( double d )  { D_D2 = d; }

  private:
    Setting *set;
    string ascii_filename;
    int ENum, ENumMax;
    int pos_type;
    double D_D2;//diameter of D2 target [mm]
    TVector3 GetUniformDistribution(double costMin=-1., double costMax = 1.);
    void ConvertToNKS2Coordinate(TVector3 &dir);
    double GetUniformR(double Rmax);
    FILE *asciiFile;
    TH1D *h_r_muon_gene, *h_theta_muon_gene;
    TH1D *h_x_muon_gene, *h_y_muon_gene, *h_z_muon_gene;
    TH2D *h2_xy_muon_gene;
    TH1D *h_cost_muon_gene;
    TH1D *h_fl_muon, *h_tof_muon;

    TVector3 pos_muon_stop, pos_at_SDD;

    //TFile *ofp;
    TLatex *latex;
    TLine *tline;
    TTree *tree_out;

    TCanvas *c[NCanvas];

};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CosmicDist::CosmicDist()
{

  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");
  gROOT->SetBatch(1);

  gStyle->SetOptDate(0);
  gStyle->SetOptFit(1);
  gStyle->SetHistFillStyle(3002);
  gStyle->SetHistFillColor(0);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetLineWidth(0);
  gStyle->SetOptDate(0);
  //gStyle->SetStatW(0.15);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatTextColor(1);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleTextColor(1);
  gStyle->SetGridWidth(0);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetLineWidth(0);
  gStyle->SetNdivisions(510); // tertiary*10000 + secondary*100 + first
  gStyle->SetOptStat("iMen");
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.13);

  //const Int_t NRGBs = 5;
  //const Int_t NCont = 255;
  //Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  //Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  //Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  //Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  //TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  //gStyle->SetNumberContours(NCont);
      
  set = new Setting();

  for(int i=0;i<NCanvas;i++){
    c[i]= new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1200,800 );
  }

  latex = new TLatex();
  tline  = new TLine();
  set->SetTLine(tline,2,2,1);

  gRandom -> SetSeed( time(NULL) ); //seed set by time


}
////////////////////////////////////////////////////////////////////////////
CosmicDist::~CosmicDist(){
  fclose(asciiFile);
  //ofp->Write();
  //ofp->Close();
}
////////////////////////////////////////////////////////////////////////////
void CosmicDist::MakeHist(string ofname){
  cout<<"MakeHist"<<endl;
  //ofp = new TFile(Form("%s",ofname.c_str()),"recreate");
  cout<<"newed ofp"<<endl;

  //ofp->cd();

  h_cost_muon_gene         = new TH1D("h_cost_muon_gene"        ,"", 200,  -1, 1);
  h_r_muon_gene            = new TH1D("h_r_muon_gene"           ,"", 100,   0, 50);
  h_theta_muon_gene        = new TH1D("h_theta_muon_gene"       ,"", 100, -1.*TMath::Pi(), TMath::Pi());
  h_x_muon_gene            = new TH1D("h_x_muon_gene"           ,"", 100, -50, 50);
  h_y_muon_gene            = new TH1D("h_y_muon_gene"           ,"", 100, -50, 50);
  h_z_muon_gene            = new TH1D("h_z_muon_gene"           ,"", 100, -50, 50);
  h_fl_muon                = new TH1D("h_fl_muon" ,"", 100, 0, 100);
  h_tof_muon               = new TH1D("h_tof_muon","", 300, 0, 300);
  h2_xy_muon_gene          = new TH2D("h2_xy_muon_gene"         ,"", 100, -50, 50, 100, -50, 50); 
  set->SetTH1(h_cost_muon_gene         ,"Cosmic cos#theta"                ,"cos#theta"    ,"Counts/0.1"    ,1,3000,0);
  set->SetTH1(h_r_muon_gene            ,"Cosmic source point (R)"         ,"R [mm]"       ,"Counts"        ,1,3000,0);
  set->SetTH1(h_theta_muon_gene        ,"Cosmic source point (#theta)"    ,"#theta [rad]" ,"Counts"        ,1,3000,0);
  set->SetTH1(h_x_muon_gene            ,"Cosmic source point X"           ,"X [mm]" ,"Counts"        ,1,3000,0);
  set->SetTH1(h_y_muon_gene            ,"Cosmic source point Y"           ,"Y [mm]" ,"Counts"        ,1,3000,0);
  set->SetTH1(h_z_muon_gene            ,"Cosmic source point Z"           ,"Z [mm]" ,"Counts"        ,1,3000,0);
  set->SetTH1(h_fl_muon                ,"Cosmic flight length"            ,"[mm]"        ,"Counts"    ,4,3000,0);
  set->SetTH1(h_tof_muon               ,"Cosmic time-of-flight"           ,"ToF [ns]"    ,"Counts"    ,4,3000,0);
  set->SetTH2(h2_xy_muon_gene          ,"Cosmic source point"             ,"X [mm]" ,"Y [mm]" );

}
////////////////////////////////////////////////////////////////////////////
void CosmicDist::Loop(){
  cout<<"Loop"<<endl;

  for(int n=0;n<ENumMax;n++){
    if(n%10000==0) cout<<n<<" / "<<ENumMax<<endl;
    double r_muBeCosmic    ;
    double theta_muBeCosmic;
    double x_muBeCosmic    ;
    double y_muBeCosmic    ;

    //////////////////////
    //hot muon simulation//
    //////////////////////
    
    //source of hot muon (Solid H2 local coordinate)
    TVector3 pos_muon, dir_muon;//muon position and direction at H2 target surface
    double r_muon     = 0.;
    double theta_muon = 0.;
    double x_muon  = 0.;
    double y_muon  = 0.;
    if(pos_type==1){
      r_muon     = GetUniformR(0.5*D_D2);
      theta_muon = 2.*TMath::Pi()*(gRandom->Uniform()-0.5);
      x_muon     = r_muon*TMath::Cos(theta_muon);
      y_muon     = r_muon*TMath::Sin(theta_muon);
    }else if(pos_type==2){
      while(1){
        x_muon     = gRandom->Gaus(0.,20.);
        y_muon     = gRandom->Gaus(0.,20.);
        r_muon     = sqrt(x_muon*x_muon + y_muon*y_muon);
        theta_muon = TMath::ATan2(x_muon, y_muon);
        if(2.*r_muon<D_D2)break;
      }
    }

    pos_muon.SetXYZ(x_muon, y_muon, 0.);
    dir_muon=GetUniformDistribution(0.);
    double cost_muon =dir_muon.Z();

    h_cost_muon_gene          ->Fill(cost_muon);
    h_r_muon_gene             ->Fill(r_muon);
    h_theta_muon_gene         ->Fill(theta_muon);
    h_x_muon_gene             ->Fill(pos_muon.X());
    h_y_muon_gene             ->Fill(pos_muon.Y());
    h_z_muon_gene             ->Fill(pos_muon.Z());
    h2_xy_muon_gene           ->Fill(pos_muon.X(), pos_muon.Y());
    cout<<Form("pos_muon.X(): %lf, pos_muon.Z(): %lf, dir_muon.X(): %lf, dir_muon.Z(): %lf, cost_muon: %lf"
          ,     pos_muon.X(),      pos_muon.Z(),      dir_muon.X(),      dir_muon.Z(),      cost_muon)<<endl;

    WriteALine(pos_muon, dir_muon, n);
  }
}
////////////////////////////////////////////////
//_______________________________________________________________________________________________________________________________

TVector3 CosmicDist::GetUniformDistribution(double costMin, double costMax)
{

  // forward peak

  //  G4double uniform_tmp; 
  //  do {
  //    uniform_tmp = G4UniformRand();
  //  } while ( G4UniformRand() >= uniform_tmp) ;
  //  G4double cost = (costMax-costMin)*uniform_tmp + costMin;

  // backford peak

  //  G4double uniform_tmp; 
  //  do {
  //    uniform_tmp = G4UniformRand();
  //  } while( G4UniformRand() <= uniform_tmp );
  //  G4double cost = (costMax-costMin)*uniform_tmp + costMin;
  

  // uniform
  double cost = (costMax-costMin)*gRandom->Uniform(1) + costMin;
  double sint = sqrt(1.-cost*cost);
  double phi  = 2.*TMath::Pi()*(gRandom->Uniform(1)-0.5);
  
  //return TVector3(cost, sint*cos(phi), sint*sin(phi));//NKS2 coordinate (x is the beam axis)
  return TVector3(sint*cos(phi), sint*sin(phi),cost);//normal coordinate (z is the beam axis)
}
//_______________________________________________________________________________________________________________________________
double CosmicDist::GetUniformR(double Rmax)
{
  double r, yy;
  
  while(1){
    r = gRandom->Uniform(Rmax);
    yy= gRandom->Uniform(1.);
    if(yy< r/Rmax)break;
  }
  
  return r;
}
//_______________________________________________________________________________________________________________________________

////////////////////////////////////////////////////////////////////////////
bool CosmicDist::Draw(){
  int ic=0;
  c[ic]->Clear();c[ic]->Divide(2,2);
  c[ic]->cd(1);gPad->SetLogz(0);  h_cost_muon_gene          ->Draw("");
  c[ic]->cd(2);gPad->SetLogz(0);  h_r_muon_gene             ->Draw("");
  c[ic]->cd(3);gPad->SetLogz(0);  h_theta_muon_gene         ->Draw("");
  c[ic]->cd(4);gPad->SetLogz(0);  h2_xy_muon_gene           ->Draw("colz");

  ic++;
  if(ic>=NCanvas)return false;
  c[ic]->Clear();c[ic]->Divide(3,2);
  c[ic]->cd(1);gPad->SetLogy(0);h2_xy_muon_gene         ->Draw("colz");
  c[ic]->cd(2);gPad->SetLogy(0);h_x_muon_gene           ->Draw("");
  c[ic]->cd(3);gPad->SetLogy(0);h_y_muon_gene           ->Draw("");
  c[ic]->cd(4);gPad->SetLogy(0);h_z_muon_gene           ->Draw("");


  return true;
}
////////////////////////////////////////////////////////////////////////////
void CosmicDist::OutputFileOpen(){
  asciiFile=fopen(ascii_filename.c_str(),"w");
  fprintf(asciiFile,"#BLTrackFile \n");
  fprintf(asciiFile,"#x y z Px Py Pz t PDGid EventID TrackID ParentID Weight\n");
  fprintf(asciiFile,"#mm mm mm MeV/c MeV/c MeV/c ns - - - - -\n");
}
////////////////////////////////////////////////////////////////////////////
void CosmicDist::WriteALine(TVector3 pos, TVector3 dir, int evID, int trID, int pareID, double weight){
  double time=0.;
  int PDGid=22;
  fprintf(asciiFile,"%.6g %.6g %.6g ",pos.X(),pos.Y(),pos.Z());
  fprintf(asciiFile,"%.6g %.6g %.6g ",dir.X(),dir.Y(),dir.Z());
  fprintf(asciiFile,"%.6g %d %d %d %d %.6g\n",time,PDGid,evID,
                                              trID,pareID,weight);
}
////////////////////////////////////////////////////////////////////////////
void CosmicDist::SaveCanvas(string PdfFileName){
  c[0]->Print(Form("%s[",PdfFileName.c_str() ));
  for(int i=0;i<NCanvas;i++)
    c[i]->Print(Form("%s" ,PdfFileName.c_str() ));
  c[NCanvas-1]->Print(Form("%s]",PdfFileName.c_str() ));
  cout<<"saved : "<<Form("%s",PdfFileName.c_str() )<<endl;
}
////////////////////////////////////////////////////////////////////////////
//////////////////////////// main //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){

  string ofname = "root/tmp.root";
  string pdfname = "pdf/hoge.pdf";
  string ofname_root;
  //string ofname_pdf = "tmp.pdf";
  bool draw_flag = true;
  double Dinput=80.;
  int ch;
  int MaxNum   = 1000;
  int DistType = 0;
  int sdd_x = 0;
  string pngname;
  extern char *optarg;
  while((ch=getopt(argc,argv,"hw:n:bp:d:t:x:"))!=-1){
    switch(ch){
    case 'w':
      ofname = optarg;
      cout<<"output filename : "<<ofname<<endl;
      break;
    case 'd':
      Dinput = atof(optarg);
      break;
    case 'n':
      MaxNum = atoi(optarg);
      break;
    case 't':
      DistType = atoi(optarg);
      break;
    case 'x':
      sdd_x = atoi(optarg);
      break;
    case 'b':
      draw_flag=false;
      cout<<"BACH MODE!"<<endl;
      break;
    case 'p':
      pdfname = optarg;
      break;
    case 'h':
      cout<<"-w : output ascii filename"<<endl;
      cout<<"-p : output pdf filename" <<endl;
      cout<<"-b : batch mode"           <<endl;
      cout<<"-d : diameter of D2 target (phi [mm])"     <<endl;
      cout<<"-t : Cosmic source poit type 0(default): from center only, 1: phi60 mm uniform, 2: gauss (sigma=20mm)"           <<endl;
      cout<<"-n : maximum number of events to be analysed "<<endl;
      return 0;
      break;
    case '?':
      cout<<"unknown option...."<<endl;
      return 0;
      break;
    default:
      cout<<"type -h to see help!!"<<endl;
      return 0;
    }
  }

  ofname_root = ofname;


  TApplication *theApp = new TApplication("App", &argc, argv);
  CosmicDist *calc = new CosmicDist();
  calc->SetFileName(ofname);
  calc->OutputFileOpen();
  calc->SetMaxEvent(MaxNum);
  calc->SetDiameter(Dinput);
  calc->SetDistType(DistType);
  calc->MakeHist("");
  //calc->SetBranch();
  calc->Loop();
  if(draw_flag){
    if(!calc->Draw())cout<<"Number of Canvas is not enough!"<<endl;
    calc->SaveCanvas(pdfname);
  }
  delete calc;
  cout<<"delete calc"<<endl;

  gSystem->Exit(1);
  theApp->Run();
  return 0;
}
