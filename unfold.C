//Funzionamento root -l unfold.C -> fa il grafico di default con c =0.25
//.x unfold.C(numero)->Fa il grafico con il valore di c specificato
//.x unfold.C(Numero_Negativo)->Salva i grafici per dei valori preimpostati: 0.3 , 0.5 , 0.7 , 0.9 , 5
const static double lambda = 420;//nm
const static double L = 20;//cm distanza sorgente
const static double a = 1.5*lambda;//dimensioni fenditura 
const static double d=6*lambda;// distanza fenditure
const static double I=1.0; //intensità sorgente

static int Eventi=5e4;
static int Bin= 100;

using namespace std;

TRandom3* r = new TRandom3();

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

string ToString(double n , int precision){
    string i=to_string((int)n);
    string d;
    double decimal=(n-(int)n)*10;
    for(int i=0; i<precision ; i++){
        d+=to_string((int)(decimal));
        decimal=(decimal-(int)decimal)*10;
    }
    return i+"."+d;
}

void unfold(double stamp = 0.25)
{
	Reset();
    gSystem->Load("/opt/RooUnfold/trunk/libRooUnfold.rootmap");
    int Cartella= system("mkdir -p UNFOLD");
	//Informazioni statistiche da stampare
	//gStyle->SetOptFit(1111);
    double x, xs;
    double c[5]={0.3 , 0.5 , 0.7 , 0.9 , 5.0};
    string c_s[5]={"0_3" , "0_5", "0_7", "0_9", "5"};
    double sigma[5]={c[0]*2*L/Bin, c[1]*2*L/Bin, c[2]*2*L/Bin, c[3]*2*L/Bin, c[4]*2*L/Bin,};
    TF1 *Diff = new TF1("Diff" , Intensita ,-L,L,0);

    int k=0 ,fine=5;

    if(stamp <= 0) {
        k=0;
        fine=5;
    }
    else{
        k=0;
        fine=1;
        sigma[0]=stamp*2*L/Bin;
        c[0]=stamp;

    }

    for(int j=k; j<fine ; j++){
        TH1D* shape = new TH1D(("Spettro c :"+ToString(c[j] , 2)).c_str() , ("Spettro c :"+ToString(c[j] , 2)).c_str() , Bin , -L , L);
        shape->SetTitleSize(50);
        shape->GetXaxis()->SetTitle("I");shape->GetXaxis()->SetTitleSize(0.055);shape->GetXaxis()->SetTitleOffset(0.7);
        shape->GetYaxis()->SetTitle("Conteggi");shape->GetYaxis()->SetTitleSize(0.055);shape->GetYaxis()->SetTitleOffset(0.85);
        // shape->SetFillColorAlpha(kYellow , 0.30);
        // shape->SetFillStyle(3008);

        TH1D* shape_smuss = new TH1D("shape_smuss" , "shape_smuss" , Bin , -L , L);
        shape_smuss->SetLineColor(kRed);
        shape_smuss->SetFillColorAlpha(kYellow , 0.30);
        shape_smuss->SetFillStyle(3004);

        TLegend* legend = new TLegend(0.1,0.7,0.4,0.9);

        RooUnfoldResponse riv (Bin , -L , L);
        for(int i=0; i<Eventi; i++){
            x=Diff->GetRandom();
            riv.Fill(smearing(x, sigma[j]),x);
        }
        
        for(int i=0; i<Eventi; i++){
            x=Diff->GetRandom();
            shape->Fill(x);
            shape_smuss->Fill(smearing(x , sigma[j]));
        }
        
        //RooUnfoldBinByBin unfold (&riv, shape_smuss);
        RooUnfoldBayes unfold (&riv, shape_smuss, 4);
        //RooUnfoldSvd unfold (&riv, shape_smuss, 21);
        TH1D* shape_reco= (TH1D*) unfold.Hreco();
        shape_reco->SetMarkerColor(kBlack);
        shape_reco->SetLineColor(kBlack);
        shape_reco->SetMarkerSize(0.7);
        shape_reco->SetMarkerStyle(20);

        int NBIN = shape_reco->GetNbinsX();
        double scarti[NBIN], errscarti[NBIN], posx[NBIN], errposx[NBIN];
        double chi=0, ev1, ev2, sigma1, sigma2;

        for(int i =0; i<NBIN; i++){
            ev1=shape->GetBinContent(i+1);
            sigma1=shape->GetBinError(i+1);

            ev2=shape_reco->GetBinContent(i+1);
            sigma2=shape_reco->GetBinError(i+1);

            //cout<<ev1<<" "<<ev2<<" "<<sigma1<<" "<<sigma2<<endl;

            if(sigma1>0 || sigma2>0) chi+=pow(ev1-ev2 , 2)/(sigma1*sigma1+sigma2*sigma2);

            scarti[i] = ev1-ev2;
            posx[i] = shape_reco->GetBinCenter(i+1);
            errscarti[i] = sqrt(sigma1*sigma1+sigma2*sigma2);
            errposx[i] = L/NBIN;
            
        }

        cout<<"Chi^2 :"<<chi<<endl;

        legend->AddEntry(shape_smuss , "Distrbuzione Osservata", "lp");
        legend->AddEntry(shape_reco , "Distribuzione ricostruita" , "lpF");
        legend->AddEntry(shape , "Distribuzione vera" ,"lp");

        gStyle->SetOptStat(0);
        TCanvas* c1 = new TCanvas();
        c1->SetGrid();
        shape->Draw();
        shape_smuss->Draw("same");
        shape_reco->Draw("Psame");
        legend->Draw();

        //TCanvas* c2 = new TCanvas();
        //c2->SetGrid();
        //unfold.Ereco().Draw("colz");

        TCanvas* c3 = new TCanvas();
        c3->SetGrid();
        TGraphErrors* gr = new TGraphErrors(NBIN,posx,scarti,errposx,errscarti);
        
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(9);
        gr->SetMarkerSize(1);
        gr->SetTitle("Unfolding Bias");
        gr->Draw("ABP");
        gr->Fit("pol0");

        
        if(stamp <=0){
            c1->SaveAs(("UNFOLD/unfold_"+c_s[j]+".png").c_str());
            c1->SaveAs(("UNFOLD/unfold_"+c_s[j]+".root").c_str());
            c1->SaveAs(("UNFOLD/unfold_"+c_s[j]+".pdf").c_str());

            //c2->SaveAs(("UNFOLD/covariance_"+c_s[j]+".png").c_str());
            //c2->SaveAs(("UNFOLD/covariance_"+c_s[j]+".root").c_str());
            //c2->SaveAs(("UNFOLD/covariance_"+c_s[j]+".pdf").c_str());

            c3->SaveAs(("UNFOLD/bias_"+c_s[j]+".png").c_str());
            c3->SaveAs(("UNFOLD/bias_"+c_s[j]+".root").c_str());
            c3->SaveAs(("UNFOLD/bias_"+c_s[j]+".pdf").c_str());

            delete c1;
            //delete c2;
            delete legend;
            delete shape_smuss;
            delete shape_reco;
            delete shape;
            riv.Delete();
            unfold.Delete();
        }
    }


    r->Delete();
    return;
}