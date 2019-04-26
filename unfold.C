
const static double lambda = 420;//nm
const static double L = 20;//cm distanza sorgente
const static double a = 1.5*lambda;//dimensioni fenditura 
const static double d=6*lambda;// distanza fenditure
const static double I=1.0; //intensità sorgente

static int Eventi=5e4;
static int Bin= 100;

using namespace std;

TRandom3* r = new TRandom3(time(0));

void Reset(){
	//Recupera la lista di tutti i canvas e li cancella
	TSeqCollection* canvases = gROOT->GetListOfCanvases();
	TIter next(canvases);
	while(TCanvas *c = (TCanvas*)next())
	{
		delete c;
	}
	//Cancella tutti gli oggetti creati
	//Necessario per evitare "warning <TROOT::Append>: Replacing existing TH1: <Name> (Potential memory leak)"
	gROOT->DeleteAll();
	return;
}

double Intensita(double* x , double *par){
    double phi = atan(x[0]/L);
    double beta= sin(phi)*d*TMath::Pi()/lambda;
    double alpha= sin(phi)*a*TMath::Pi()/lambda;
    double i=4*I*pow(sin(alpha)*cos(beta)/alpha, 2);
    return i;
}

double smearing(double x, double sigma){
    return x+(r->Gaus(0 , sigma));
}

void unfold()
{
	Reset();
	//Informazioni statistiche da stampare
	//gStyle->SetOptFit(1111);
    double x, xs,sigma;
    double c=5;
    sigma=c*2*L/Bin;
    TF1 *Diff = new TF1("Diff" , Intensita ,-L,L,0);

    TH1D* shape = new TH1D("shape" , "shape" , Bin , -L , L);
    // shape->SetFillColorAlpha(kYellow , 0.30);
    // shape->SetFillStyle(3008);

    TH1D* shape_smuss = new TH1D("shape_smuss" , "shape_smuss" , Bin , -L , L);
    shape_smuss->SetLineColor(kRed);
    shape_smuss->SetFillColorAlpha(kYellow , 0.30);
    shape_smuss->SetFillStyle(3004);

    gSystem->Load("/opt/RooUnfold/trunk/libRooUnfold.rootmap");

    RooUnfoldResponse riv (Bin , -L , L);

    for(int i=0; i<Eventi; i++){
        x=Diff->GetRandom();
        riv.Fill(smearing(x, sigma),x);
    }
    
    for(int i=0; i<Eventi; i++){
        x=Diff->GetRandom();
        shape->Fill(x);
        shape_smuss->Fill(smearing(x , sigma));
    }
    RooUnfoldBinByBin unfold (&riv, shape_smuss);
    //RooUnfoldBayes deconvoluzione (&riv, shape_smuss, 4);
    TH1D* shape_reco= (TH1D*) unfold.Hreco();
    shape_reco->SetMarkerColor(kBlack);
    shape_reco->SetMarkerSize(0.7);
    shape_reco->SetMarkerStyle(20);

    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas();
    c1->SetGrid();
    shape->Draw();
    shape_smuss->Draw("same");
    shape_reco->Draw("Psame");
    
    


    r->Delete();
    riv.Delete();
    unfold.Delete();
    return;
}