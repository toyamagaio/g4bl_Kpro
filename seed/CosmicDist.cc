#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <cstdlib>
#include <cerrno>
#include <cmath>
#include <getopt.h>
#include <algorithm>
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
  void SetAnglePower( double n ) { angle_power = n; }
  void SetMomentumDistType( int type ) { momentum_type = type; }
  void SetMomentumRange( double min, double max ) { momentum_min = min; momentum_max = max; }
  void SetMomentumCoefficients( double c0, double c1, double c2 ) {
    momentum_c0 = c0; momentum_c1 = c1; momentum_c2 = c2;
  }
  void SetPositionShape( int shape ) { position_shape = shape; }
  void SetPositionSize( double width, double height ) {
    position_width = width; position_height = height;
  }
  bool ValidateConfiguration( string &error ) const;

  private:
    Setting *set;
    string ascii_filename;
    int ENum, ENumMax;
    int pos_type;
    double D_D2;//diameter of D2 target [mm]
    double angle_power;
    int momentum_type; // 0: uniform, 1: linear, 2: quadratic
    double momentum_min, momentum_max;
    double momentum_c0, momentum_c1, momentum_c2;
    int position_shape; // 0: circle, 1: rectangle
    double position_width, position_height;
    TVector3 GetCosPowerDistribution(double n);
    double GetMomentum();
    double MomentumShape(double x) const;
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
  : ENum(0), ENumMax(0), pos_type(0), D_D2(80.), angle_power(0.),
    momentum_type(0), momentum_min(1.), momentum_max(1.),
    momentum_c0(1.), momentum_c1(0.), momentum_c2(0.),
    position_shape(0), position_width(80.), position_height(80.), asciiFile(0)
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
  if(asciiFile) fclose(asciiFile);
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
      if(position_shape==1){
        x_muon = gRandom->Uniform(-0.5*position_width, 0.5*position_width);
        y_muon = gRandom->Uniform(-0.5*position_height, 0.5*position_height);
        r_muon = sqrt(x_muon*x_muon + y_muon*y_muon);
        theta_muon = TMath::ATan2(y_muon, x_muon);
      }else{
        r_muon     = GetUniformR(0.5*D_D2);
        theta_muon = 2.*TMath::Pi()*(gRandom->Uniform()-0.5);
        x_muon     = r_muon*TMath::Cos(theta_muon);
        y_muon     = r_muon*TMath::Sin(theta_muon);
      }
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
    dir_muon=GetCosPowerDistribution(angle_power);
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

    WriteALine(pos_muon, dir_muon*GetMomentum(), n);
  }
}
////////////////////////////////////////////////
//_______________________________________________________________________________________________________________________________

TVector3 CosmicDist::GetCosPowerDistribution(double n)
{
  // A flux proportional to cos^n(theta) per unit solid angle has
  // p(cos(theta)) proportional to cos^n(theta).  Invert its CDF on
  // the downward hemisphere: cos(theta) = U^(1/(n+1)).
  double cost = pow(gRandom->Uniform(), 1.0/(n+1.0));
  double sint = sqrt(1.-cost*cost);
  double phi  = 2.*TMath::Pi()*(gRandom->Uniform()-0.5);
  
  //return TVector3(cost, sint*cos(phi), sint*sin(phi));//NKS2 coordinate (x is the beam axis)
  return TVector3(sint*cos(phi), sint*sin(phi),cost);//normal coordinate (z is the beam axis)
}
//_______________________________________________________________________________________________________________________________
double CosmicDist::MomentumShape(double x) const
{
  double value = momentum_c0 + momentum_c1*x;
  if(momentum_type==2) value += momentum_c2*x*x;
  return value;
}
//_______________________________________________________________________________________________________________________________
double CosmicDist::GetMomentum()
{
  if(momentum_min==momentum_max)
    return momentum_min;
  if(momentum_type==0)
    return gRandom->Uniform(momentum_min,momentum_max);

  // Invert the normalized polynomial CDF on x in [0,1]. Bisection avoids
  // special-case roots and remains stable for linear and quadratic shapes.
  double norm = momentum_c0 + 0.5*momentum_c1;
  if(momentum_type==2) norm += momentum_c2/3.0;
  double target = gRandom->Uniform()*norm;
  double low = 0.0, high = 1.0;
  for(int i=0;i<48;i++){
    double x = 0.5*(low+high);
    double integral = momentum_c0*x + 0.5*momentum_c1*x*x;
    if(momentum_type==2) integral += (momentum_c2/3.0)*x*x*x;
    if(integral < target) low = x;
    else high = x;
  }
  return momentum_min + (momentum_max-momentum_min)*0.5*(low+high);
}
//_______________________________________________________________________________________________________________________________
bool CosmicDist::ValidateConfiguration(string &error) const
{
  if(!std::isfinite(angle_power) || angle_power<0.){
    error = "angle power n must be finite and non-negative";
    return false;
  }
  if(ENumMax<0){
    error = "number of events must be non-negative";
    return false;
  }
  if(momentum_type<0 || momentum_type>2){
    error = "momentum distribution must be uniform, linear, or quadratic";
    return false;
  }
  if(!std::isfinite(momentum_min) || !std::isfinite(momentum_max) ||
     momentum_min<0. || momentum_max<momentum_min){
    error = "momentum range must satisfy 0 <= min <= max";
    return false;
  }
  if(momentum_type!=0 && momentum_min!=momentum_max){
    if(!std::isfinite(momentum_c0) || !std::isfinite(momentum_c1) ||
       (momentum_type==2 && !std::isfinite(momentum_c2))){
      error = "momentum coefficients must be finite";
      return false;
    }
    double min_value = MomentumShape(0.);
    min_value = std::min(min_value, MomentumShape(1.));
    if(momentum_type==2 && momentum_c2>0.){
      double vertex = -momentum_c1/(2.*momentum_c2);
      if(vertex>0. && vertex<1.) min_value = std::min(min_value, MomentumShape(vertex));
    }
    double norm = momentum_c0 + 0.5*momentum_c1;
    if(momentum_type==2) norm += momentum_c2/3.0;
    if(min_value<0. || !std::isfinite(norm) || norm<=0.){
      error = "momentum polynomial must be non-negative on [0,1] and have positive integral";
      return false;
    }
  }
  if(position_shape<0 || position_shape>1){
    error = "position shape must be circle or rectangle";
    return false;
  }
  if(pos_type<0 || pos_type>2){
    error = "position type must be 0, 1, or 2";
    return false;
  }
  if(pos_type==1 && position_shape==1 &&
     (!std::isfinite(position_width) || !std::isfinite(position_height) ||
      position_width<=0. || position_height<=0.)){
    error = "rectangle width and height must be positive";
    return false;
  }
  if((pos_type==1 && position_shape==0) || pos_type==2){
    if(!std::isfinite(D_D2) || D_D2<=0.){
      error = "circular source diameter must be positive";
      return false;
    }
  }
  return true;
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
bool ParseDoubleOption(const char *text, double &value)
{
  char *end = 0;
  errno = 0;
  value = strtod(text, &end);
  return end!=text && *end=='\0' && errno!=ERANGE && std::isfinite(value);
}

bool ParseIntOption(const char *text, int &value)
{
  char *end = 0;
  errno = 0;
  long parsed = strtol(text, &end, 10);
  if(end==text || *end!='\0' || errno==ERANGE || parsed<(-2147483647L-1) || parsed>2147483647L)
    return false;
  value = static_cast<int>(parsed);
  return true;
}

bool ParseMomentumType(const char *text, int &type)
{
  string value(text);
  if(value=="uniform" || value=="0") type=0;
  else if(value=="linear" || value=="1") type=1;
  else if(value=="quadratic" || value=="2") type=2;
  else return false;
  return true;
}

bool ParsePositionShape(const char *text, int &shape)
{
  string value(text);
  if(value=="circle" || value=="circular") shape=0;
  else if(value=="rectangle" || value=="rectangular" || value=="rect") shape=1;
  else return false;
  return true;
}

int main(int argc, char** argv){

  string ofname = "root/tmp.root";
  string pdfname = "pdf/hoge.pdf";
  string ofname_root;
  //string ofname_pdf = "tmp.pdf";
  bool draw_flag = true;
  double Dinput=80.;
  double anglePower=0.;
  int MomentumType=0;
  double MomentumMin=1., MomentumMax=1.;
  double MomentumC0=1., MomentumC1=0., MomentumC2=0.;
  int PositionShape=0;
  double PositionWidth=80., PositionHeight=80.;
  int ch;
  int MaxNum   = 1000;
  int DistType = 0;
  int sdd_x = 0;
  string pngname;
  extern char *optarg;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"output", required_argument, 0, 'w'},
    {"pdf", required_argument, 0, 'p'},
    {"batch", no_argument, 0, 'b'},
    {"diameter", required_argument, 0, 'd'},
    {"position-type", required_argument, 0, 't'},
    {"events", required_argument, 0, 'n'},
    {"angle-power", required_argument, 0, 'a'},
    {"momentum", required_argument, 0, 'm'},
    {"momentum-min", required_argument, 0, 'l'},
    {"momentum-max", required_argument, 0, 'u'},
    {"p0", required_argument, 0, 'c'},
    {"p1", required_argument, 0, 'A'},
    {"p2", required_argument, 0, 'B'},
    {"position-shape", required_argument, 0, 'g'},
    {"width", required_argument, 0, 'W'},
    {"height", required_argument, 0, 'H'},
    {0, 0, 0, 0}
  };
  while((ch=getopt_long(argc,argv,"hw:n:bp:d:t:x:a:m:l:u:c:A:B:g:W:H:",long_options,0))!=-1){
    switch(ch){
    case 'w':
      ofname = optarg;
      cout<<"output filename : "<<ofname<<endl;
      break;
    case 'd':
      if(!ParseDoubleOption(optarg,Dinput)){ cerr<<"invalid diameter: "<<optarg<<endl; return 1; }
      break;
    case 'n':
      if(!ParseIntOption(optarg,MaxNum)){ cerr<<"invalid event count: "<<optarg<<endl; return 1; }
      break;
    case 't':
      if(!ParseIntOption(optarg,DistType)){ cerr<<"invalid position type: "<<optarg<<endl; return 1; }
      break;
    case 'x':
      if(!ParseIntOption(optarg,sdd_x)){ cerr<<"invalid x value: "<<optarg<<endl; return 1; }
      break;
    case 'a':
      if(!ParseDoubleOption(optarg,anglePower)){ cerr<<"invalid angle power: "<<optarg<<endl; return 1; }
      break;
    case 'm':
      if(!ParseMomentumType(optarg,MomentumType)){ cerr<<"momentum type must be uniform, linear, or quadratic"<<endl; return 1; }
      break;
    case 'l':
      if(!ParseDoubleOption(optarg,MomentumMin)){ cerr<<"invalid minimum momentum: "<<optarg<<endl; return 1; }
      break;
    case 'u':
      if(!ParseDoubleOption(optarg,MomentumMax)){ cerr<<"invalid maximum momentum: "<<optarg<<endl; return 1; }
      break;
    case 'c':
      if(!ParseDoubleOption(optarg,MomentumC0)){ cerr<<"invalid momentum coefficient p0: "<<optarg<<endl; return 1; }
      break;
    case 'A':
      if(!ParseDoubleOption(optarg,MomentumC1)){ cerr<<"invalid momentum coefficient p1: "<<optarg<<endl; return 1; }
      break;
    case 'B':
      if(!ParseDoubleOption(optarg,MomentumC2)){ cerr<<"invalid momentum coefficient p2: "<<optarg<<endl; return 1; }
      break;
    case 'g':
      if(!ParsePositionShape(optarg,PositionShape)){ cerr<<"position shape must be circle or rectangle"<<endl; return 1; }
      break;
    case 'W':
      if(!ParseDoubleOption(optarg,PositionWidth)){ cerr<<"invalid rectangle width: "<<optarg<<endl; return 1; }
      break;
    case 'H':
      if(!ParseDoubleOption(optarg,PositionHeight)){ cerr<<"invalid rectangle height: "<<optarg<<endl; return 1; }
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
      cout<<"-t : source position type 0(default): center, 1: uniform area, 2: truncated Gaussian (sigma=20 mm)"<<endl;
      cout<<"-n : maximum number of events to be analysed "<<endl;
      cout<<"-a, --angle-power : flux per solid angle proportional to cos(theta)^n on downward hemisphere (default n=0)"<<endl;
      cout<<"-m, --momentum : momentum distribution: uniform, linear, or quadratic (default uniform)"<<endl;
      cout<<"-l, --momentum-min : minimum momentum in MeV/c (default 1)"<<endl;
      cout<<"-u, --momentum-max : maximum momentum in MeV/c (default 1)"<<endl;
      cout<<"-c, --p0 : constant momentum-polynomial coefficient (default 1)"<<endl;
      cout<<"-A, --p1 : linear momentum-polynomial coefficient (default 0)"<<endl;
      cout<<"-B, --p2 : quadratic momentum-polynomial coefficient (default 0)"<<endl;
      cout<<"             f(x)=p0+p1*x (+p2*x^2), x=(p-pmin)/(pmax-pmin); density is normalized automatically"<<endl;
      cout<<"-g, --position-shape : circle or rectangle for uniform-area position type -t 1 (default circle)"<<endl;
      cout<<"-W, --width : rectangle width in mm (default 80)"<<endl;
      cout<<"-H, --height : rectangle height in mm (default 80)"<<endl;
      cout<<"Example: -a 2 -m quadratic -l 100 -u 1000 -c 1 -A 0 -B 2"<<endl;
      cout<<"          -t 1 -g rectangle -W 120 -H 80"<<endl;
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

  TApplication *theApp = new TApplication("App", &argc, argv);
  string configError;
  CosmicDist *calc = new CosmicDist();
  calc->SetMaxEvent(MaxNum);
  calc->SetDiameter(Dinput);
  calc->SetDistType(DistType);
  calc->SetAnglePower(anglePower);
  calc->SetMomentumDistType(MomentumType);
  calc->SetMomentumRange(MomentumMin,MomentumMax);
  calc->SetMomentumCoefficients(MomentumC0,MomentumC1,MomentumC2);
  calc->SetPositionShape(PositionShape);
  calc->SetPositionSize(PositionWidth,PositionHeight);
  if(!calc->ValidateConfiguration(configError)){
    cerr<<"invalid configuration: "<<configError<<endl;
    delete calc;
    return 1;
  }

  ofname_root = ofname;


  calc->SetFileName(ofname);
  calc->OutputFileOpen();
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
