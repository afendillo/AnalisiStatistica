
R__LOAD_LIBRARY(string)

static const double Pi= TMath::Pi();
static const int NFibre= 32;
static const double spessore= 0.3; //cm
static const double lunghezza= 10; //cm
static const double FAC= 0.04;
static const double Q_e= 1.60217662e-19;
static const double L= 210; //cm
static const double I0= 1;//cm-2s-1
static const double Conv= 180/Pi;
static const double EGamma=197.326e-6*2*Pi/420; //mev nm-->h tagliato
static double G=1e6;//Gain
static double QE=0.5;//efficienza di conversione


using namespace std;

struct posizione3D{double x , y , theta, phi;int piano, fibra, flag=0;};

struct segnale {vector <double> Nf, q, R; int flag;};

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

void drawtext()
{
   Int_t i,n;
   Double_t x,y;
   TLatex l;
   l.SetTextSize(0.025);
   l.SetTextFont(42);
   l.SetTextAlign(21);
   l.SetTextColor(kBlue);
   TGraph *g = (TGraph*)gPad->GetListOfPrimitives()->FindObject("Graph");
   n = g->GetN();
   for (i=0; i<n; i++) {
      g->GetPoint(i,x,y);
      l.PaintText(x,y+0.2,Form("(%4.2f,%4.2f)",x,y));
   }
}


double FibraPos(posizione3D pos){
    double fibra;
    if (pos.theta!=Pi/2)
    {
        if((pos.x/spessore)-(int)(pos.x/spessore)==0 && pos.phi>Pi )
        {
            fibra=floor(pos.x/spessore);
        }
        else fibra=ceil(pos.x/spessore);//fibra colpita
    }

    if (fibra>NFibre || fibra <=0) fibra = NFibre*(fibra>NFibre)+(fibra<=0);

    return fibra;

}
double* Angle(TRandom3* rdm)
{
    double * angle= new double[2];
    angle[0]=(rdm->Rndm()+1)*Pi/2;//theta
    double hit=rdm->Rndm();
    while(hit>pow(cos(angle[0]), 2)){
        hit=rdm->Rndm();
        angle[0]=(rdm->Rndm()+1)*Pi/2;
    }
    angle[1]=rdm->Rndm()*2*Pi;//phi
    
    return angle;
}

posizione3D CoordPoint1P(TRandom3* rdm)
{
    //origine estremo del piano quadrato-->prima fibra = pos 0.3
    posizione3D pos;
    pos.x=rdm->Rndm()*(NFibre)*spessore;
    pos.y=rdm->Rndm()*lunghezza;
    double *posTP=Angle(rdm);
    pos.theta=posTP[0];
    pos.phi=posTP[1];
    pos.fibra=FibraPos(pos);
    delete posTP;
    
    return pos;
}

posizione3D InversionXY(posizione3D pos){
    posizione3D NewPos=pos;
    NewPos.fibra=FibraPos(NewPos);
    if(pos.piano==2){
        NewPos.x=pos.y;//inverto per poter calcolare correttamente la posizione della fibra
        NewPos.fibra=FibraPos(NewPos);
        NewPos.x=pos.x;//Li riporto a normalità per poter sfruttare lo stesso calcolo nella proiezione
    }
    return NewPos;
}


posizione3D Proiezione(posizione3D pos)
{
    posizione3D NewPos;
    NewPos=pos;
    if (pos.theta!=Pi/2)
    {
        NewPos.x=pos.x-1*sin(NewPos.phi)*spessore*tan(NewPos.theta);//Asse Z positivo verso l'alto, i muoni vanno verso il basso
        NewPos.y=pos.y-1*cos(NewPos.phi)*spessore*tan(NewPos.theta);
        NewPos.fibra=(InversionXY(NewPos)).fibra;
        //scelta coord x=R*sin(phi)*sin(theta), y=R*cos(phi)*sin(theta), z=R*cos(theta)--->le coord vanno messe con segno
    }

    if(NewPos.x>10 || NewPos.x<0 || NewPos.y<0 || NewPos.y>9.6) NewPos.flag=1;
    return NewPos;
}


double QRaccolta(double x, double E){
    //x distanza pmt, E energia depositata
    double nGamma=E/EGamma;
    return nGamma*FAC*exp(-x/L)*Q_e*G*QE;
}

double Uniforme (double x1, double x2){

    TRandom3* rand=new TRandom3(time(0));
    double max, min , result;
    max=std::max(x1,x2);
    min=std::min(x1,x2);
    result=rand->Rndm()*(max-min)+min;
    rand->Delete();
    return result;
}

segnale deposito1P(posizione3D pos1 , posizione3D pos2, TH1D* hist){

    segnale sig;
    //TGraph2D* graph = new TGraph2D();
    double phi= pos1.phi, yout, dist , y=pos1.y , x=pos1.x, z=0, NewX, zout, Carica;
    int NumFibre= abs(pos1.fibra-pos2.fibra), k,j=0;
    k=1*(pos2.x>pos1.x)-1*(pos2.x<pos1.x);//per sapere se andare a destra (x cresce) o sinistra (z decresce)
    int corr=-1*(pos2.x<pos1.x);

    if(pos1.piano*pos2.piano==1){//per il piano 1
            // graph->SetPoint(j, x , y , z);
            // cout<<"----------------------PIANO 1------------------------------\n";
            // cout<<"X: "<<x<<" Y: " <<y <<" Z: "<<z <<endl;
            for(int i=pos1.fibra+corr; x<spessore*NFibre && y<lunghezza; i=i+k){
                //Calcola la nuova X=edge della fibra o x finale
                //NOTA: se la x finale è fuori accettanza verrà usata l'edge della fibra finale
                if ((i-corr)!=pos2.fibra) NewX=i*spessore;
                else if((i-corr)==pos2.fibra) NewX=pos2.x*(pos2.x<9.6)+i*spessore*(pos2.x>9.6);

                yout=y+(NewX-x)/tan(phi);//nuova Y
                zout=z+(NewX-x)/(sin(phi)*tan(pos1.theta));

                dist=sqrt(pow(NewX-x, 2)+pow(yout-y,2)+pow(zout-z,2));//dist per deposito energia

                //salvo i vari dati: energie depositata, carica, fibra e piano
                //(sig.en).push_back(dist*1.89/EGamma);
                Carica=QRaccolta(Uniforme(y,yout),dist*1.89);
                (sig.q).push_back(Carica);//1.89MeV/cm densità energia media depositata
                (sig.Nf).push_back(i-corr);//le fibre sono definite in base alla loro UpperEdge
                sig.flag=pos1.piano;
                (sig.R).push_back(sqrt(pow(NewX-pos1.x, 2)+pow(yout-pos2.y,2)+pow(zout,2)));

                //aggiorno le coordinate per il calcolo del punto successivo
                y=yout;
                z=zout;
                x=NewX;

                // cout<<"X: "<<x<<" Y: " <<y <<" Z: "<<z <<" R: "<<dist<<endl;
                // j++;
                // graph->SetPoint(j, x , y , z);
                hist->Fill(i-corr , Carica);
                if((i-corr)==pos2.fibra) break;
                
            }
    }

    // TCanvas* c = new TCanvas();

    // graph->SetTitle("Spost 1P; X; Y; Z ");
    // graph->SetMarkerStyle(20);
    // graph->Draw("LINE PL");
    return sig;
}


segnale deposito2P(posizione3D pos1 , posizione3D pos2 , TH1D* hist){

    segnale sig;
    //TGraph2D* graph = new TGraph2D();
    double phi= pos1.phi, xout, dist , y=pos1.y , x=pos1.x, z=-0.3, NewY, zout , Carica;
    int NumFibre= abs(pos1.fibra-pos2.fibra), k,j=0;
    k=1*(pos2.y>pos1.y)-1*(pos2.y<pos1.y);//per sapere se andare a destra (x cresce) o sinistra (z decresce)
    int corr=-1*(pos2.y<pos1.y);

    if(pos1.piano*pos2.piano==4){//per il piano 2
            // graph->SetPoint(j, x , y , z);
            // cout<<"----------------------PIANO 2------------------------------\n";
            // cout<<"X: "<<x<<" Y: " <<y <<" Z: "<<z <<endl;
            for(int i=pos1.fibra+corr; y<spessore*NFibre && x<lunghezza; i=i+k){
                //Calcola la nuova Y=edge della fibra o y finale
                //NOTA: se la y finale è fuori accettanza verrà usata l'edge della fibra finale
                if ((i-corr)!=pos2.fibra) NewY=i*spessore;
                else if((i-corr)==pos2.fibra) NewY=pos2.y*(pos2.y<9.6)+i*spessore*(pos2.y>9.6);

                xout=x+(NewY-y)*tan(phi);//nuova x
                zout=z+(NewY-y)/(cos(phi)*tan(pos1.theta));//nuovo z

                dist=sqrt(pow(xout-x, 2)+pow(NewY-y,2)+pow(zout-z,2));//dist per deposito energia

                //salvo i vari dati: energie depositata, carica, fibra e piano
                //(sig.en).push_back(dist*1.89/EGamma);
                Carica=QRaccolta(Uniforme(x,xout),dist*1.89);
                (sig.q).push_back(Carica);//1.89MeV/cm densità energia media depositata
                (sig.Nf).push_back(i-corr);//le fibre sono definite in base alla loro UpperEdge
                sig.flag=pos1.piano;
                (sig.R).push_back(sqrt(pow(xout-pos1.x, 2)+pow(NewY-pos1.y,2)+pow(zout+0.3,2)));

                //aggiorno le coordinate per il calcolo del punto successivo
                y=NewY;
                z=zout;
                x=xout;

                // cout<<"X: "<<x<<" Y: " <<y <<" Z: "<<z <<" R: "<<dist<<endl;
                // j++;
                // graph->SetPoint(j, x , y , z);
                hist->Fill(i-corr, Carica);
                if((i-corr)==pos2.fibra) break;
                
            }
    }
    // TCanvas* c = new TCanvas();

    // graph->SetTitle("Spost 2P; X; Y; Z ");
    // graph->SetMarkerStyle(20);
    // graph->Draw("LINE PL");
    return sig;
}

void StampaSignal(segnale sig){
    for(int i =0 ; i<(sig.q).size(); i++){
        cout<<"Piano: "<<sig.flag<<" , Fibra: "<<sig.Nf[i]<<" , Carica: "<<sig.q[i]<<" , R: "<<sig.R[i]<<endl;
        cout<<"-----------------------------------------------------------\n";
    }
    return;
}

posizione3D SetCoord(double x , double y, double theta, double phi){
    posizione3D pos;
    pos.x=x;
    pos.y=y;
    pos.theta=theta/Conv;//Theta e phi angoli in gradi
    pos.phi=phi/Conv;
    pos.fibra=FibraPos(pos);
    return pos;
}

void FAMU(int stamp=0){
    Reset();

    int NMuon=1e6;

    TRandom3* rand= new TRandom3();
    posizione3D posPiano1Up, posPiano1D , posPiano2D, posPiano2Up;
    segnale sig, sig2;

    TH1D* Htheta = new TH1D("Theta" , "Theta", 100  , Pi/2 ,Pi);
    Htheta->GetXaxis()->SetTitle("#phi");Htheta->GetXaxis()->SetTitleSize(0.055);Htheta->GetXaxis()->SetTitleOffset(0.7);
    // Htheta->GetYaxis()->SetTitle("Conteggi");Htheta->GetYaxis()->SetTitleSize(0.055);Htheta->GetYaxis()->SetTitleOffset(0.85);

    double Gain[3] = {5*1e6 , 1e6 , 5*1e5};
    double Eff[2] = {0.5 , 0.3};
    string nameG[3]={"G: 5x10^{6}", "G: 10^{6}","G 5x10^{5}"};
    string nameE[2]={" Q_{E}: 50%", " Q_{E}: 30%"};

    TH1D* h1[6];
    TH1D* h2[6];

    
    TH1D* HPhi = new TH1D("Phi" , "Phi", 100  , 0 ,2*Pi);
    HPhi->GetXaxis()->SetTitle("#phi");HPhi->GetXaxis()->SetTitleSize(0.055);HPhi->GetXaxis()->SetTitleOffset(0.7);
    // HPhi->GetYaxis()->SetTitle("Conteggi");HPhi->GetYaxis()->SetTitleSize(0.055);HPhi->GetYaxis()->SetTitleOffset(0.85);

    TH2D* HPos1 = new TH2D("Pos1" , "Pos1", NFibre  , 0,NFibre*spessore, NFibre ,0 ,lunghezza);
    HPos1->SetTitle("Distribuzione Piano 1 (x,y); X; Y;  ");HPos1->SetTitleSize(50);
    HPos1->GetXaxis()->SetTitleSize(0.055);HPos1->GetXaxis()->SetTitleOffset(0.7);
    HPos1->GetYaxis()->SetTitleSize(0.055);HPos1->GetYaxis()->SetTitleOffset(0.7);

    TH2D* HPos2 = new TH2D("Pos2" , "Pos2", NFibre  , 0, lunghezza, NFibre ,0 , NFibre*spessore);
    HPos2->SetTitle("Distribuzione Piano 2 (x,y); X; Y;  ");HPos2->SetTitleSize(50);
    HPos2->GetXaxis()->SetTitleSize(0.055);HPos2->GetXaxis()->SetTitleOffset(0.7);
    HPos2->GetYaxis()->SetTitleSize(0.055);HPos2->GetYaxis()->SetTitleOffset(0.7);

    double xl=0.65 , yl=0.13;
	
	auto legend = new TLegend(xl,yl,xl+0.25 ,yl+0.25);
    legend->SetHeader("Gain e Efficienza" , "c");

    auto legend2 = new TLegend(xl,yl,xl+0.25 ,yl+0.25);
    legend2->SetHeader("Gain e Efficienza" , "c");

    int hc=-1;
    for(int j=0 ; j<3;j++){
        G=Gain[j];
        for(int k=0 ; k<2;k++){
            hc++;
            h1[hc]=new TH1D( ("Primo Piano " +nameG[j]+nameE[k]).c_str(),("Primo Piano " +nameG[j]+nameE[k]).c_str(), NFibre  ,0.5 ,NFibre+0.5);
            h1[hc]->SetTitleSize(50);
            h1[hc]->GetXaxis()->SetTitle("# fibra");h1[hc]->GetXaxis()->SetTitleSize(0.055);h1[hc]->GetXaxis()->SetTitleOffset(0.7);
            h1[hc]->GetYaxis()->SetTitle("carica");h1[hc]->GetYaxis()->SetTitleSize(0.055);h1[hc]->GetYaxis()->SetTitleOffset(1.1);
            h1[hc]->SetLineColor(j+2);
            h1[hc]->SetLineStyle(k+1);

            h2[hc]=new TH1D( ("Secondo Piano " +nameG[j]+nameE[k]).c_str(),("Secondo Piano " +nameG[j]+nameE[k]).c_str(), NFibre  ,0.5 ,NFibre+0.5);
            h2[hc]->SetTitleSize(50);
            h2[hc]->GetXaxis()->SetTitle("# fibra");h2[hc]->GetXaxis()->SetTitleSize(0.055);h2[hc]->GetXaxis()->SetTitleOffset(0.7);
            h2[hc]->GetYaxis()->SetTitle("carica");h2[hc]->GetYaxis()->SetTitleSize(0.055);h2[hc]->GetYaxis()->SetTitleOffset(1.1);
            h2[hc]->SetLineColor(j+2);
            h2[hc]->SetLineStyle(k+1);

            QE=Eff[k];
            cout<<"Ciclo: "<<hc+1<<endl;
            for(int i=0; i<NMuon;i++){
                posPiano1Up = CoordPoint1P(rand);
                //posPiano1Up=SetCoord(0.256784 , 0.32698 , 105.789 , 45.105);//setta le coord a mano, angoli in gradi (90<theta<180, phi<360)
                posPiano1Up.piano=1;
                
                posPiano1D= Proiezione(posPiano1Up);

                if(i%(int)(NMuon/100)==0) 
                {
                    for (int k = 0 ; k<i/(int)(NMuon/100)+1; k++) cout << "*" ;
                    cout<<" "<< i/(int)(NMuon/100)+1<< "% \r";
                    cout<<flush;
                }


                if (posPiano1D.flag==0) {
                    sig = deposito1P(posPiano1Up, posPiano1D , h1[hc]);//salva in sig la Carica raccolta, numero fibra, piano e R.
                    posPiano2Up=posPiano1D;
                    posPiano2Up.piano=2;
                    posPiano2Up =InversionXY(posPiano2Up);
                    
                    posPiano2D=Proiezione(posPiano2Up);
                    posPiano2D.piano=2;

                    sig2=deposito2P(posPiano2Up,posPiano2D , h2[hc]);
                    if(hc==1){
                        HPos1->Fill(posPiano1Up.x , posPiano1Up.y);//Punti d'entrata-->uniforme.
                        HPos2->Fill(posPiano2Up.x , posPiano2Up.y);
                    }
                }
                if(hc==1){
                    Htheta->Fill(posPiano1Up.theta);
                    HPhi->Fill(posPiano1Up.phi);
                }
                
            }
            cout<<"\n";
	        legend->AddEntry(h1[hc],(nameG[j]+nameE[k]).c_str(),"lp");
	        legend2->AddEntry(h2[hc],(nameG[j]+nameE[k]).c_str(),"lp");
        }
    }
    cout<<endl;

    double m1=0,M1=0, m2=0,M2=0;
    m1=h1[0]->GetBinContent(h1[0]->GetMinimumBin());
    m2=h2[0]->GetBinContent(h2[0]->GetMinimumBin());
    for (int i=0; i<6;i++){
        if(M1<h1[i]->GetBinContent(h1[i]->GetMaximumBin())) M1=h1[i]->GetBinContent(h1[i]->GetMaximumBin());
        if(m1>h1[i]->GetBinContent(h1[i]->GetMinimumBin())) m1=h1[i]->GetBinContent(h1[i]->GetMinimumBin());

        if(M2<h2[i]->GetBinContent(h1[i]->GetMaximumBin())) M2=h2[i]->GetBinContent(h2[i]->GetMaximumBin());
        if(m2>h2[i]->GetBinContent(h2[i]->GetMinimumBin())) m2=h2[i]->GetBinContent(h2[i]->GetMinimumBin());
    }
    rand->Delete();

    gStyle->SetOptStat(0);

    TCanvas* c = new TCanvas("Distribuzione #theta" , "Distribuzione #theta");
    c->SetGrid();
    c->SetCanvasSize(600, 600);
    c->SetWindowSize(600+8, 600+28);
    //c->SetWindowSize(1500 , 780);
    Htheta->Draw();

    int wh=600 , h=600;

    TCanvas* c1 = new TCanvas("Distribuzione #phi", "Distribuzione #phi");
    c1->SetGrid();
    c1->SetCanvasSize(600, 600);
    c1->SetWindowSize(600+8, 600+28);
    //c1->SetWindowSize(1500 , 780);
    HPhi->GetYaxis()->SetRangeUser(0 , HPhi->GetBinContent(HPhi->GetMaximumBin())*1.2);
    HPhi->Draw();

    TCanvas* c2 = new TCanvas("Distribuzione Posizione piano 1" , "Distribuzione Posizione piano 1");
    //c2->SetGrid();
    c2->SetCanvasSize(wh, h);
    c2->SetWindowSize(wh+28, h+28);
    c2->SetRightMargin(0.13);
    HPos1->Draw("COLZ");

    TCanvas* c3 = new TCanvas("Distribuzione Posizione piano 2" , "Distribuzione Posizione piano 2");
    //c3->SetGrid();
    c3->SetCanvasSize( wh, h);
    c3->SetWindowSize(wh+28, h+28);
    c3->SetRightMargin(0.13);
    HPos2->Draw("COLZ");

    h1[0]->SetTitle("Primo piano");
    h1[0]->SetName("Primo piano");

    h2[0]->SetTitle("Secondo piano");
    h2[0]->SetName("Primo piano");

    TCanvas* c4= new TCanvas("Primo Piano" , "Primo Piano");
    c4->SetGrid();
    c4->SetCanvasSize(wh, h);
    c4->SetWindowSize(wh+28, h+28);
    c4->SetLeftMargin(0.13);
    c4-> SetLogy();
    h1[0]->GetYaxis()->SetRangeUser(m1/10 , M1*2);
    //h1[0]->SetStats(0);
    for(int i=0; i<6; i++) h1[i]->Draw("HIST same");
    legend->Draw();

    TCanvas* c5 = new TCanvas("Secondo Piano" , "Secondo Piano");
    c5->SetGrid();
    c5-> SetLogy();
    c5->SetCanvasSize( wh, h);
    c5->SetWindowSize(wh+28, h+28);
    c5->SetLeftMargin(0.13);
    
    h2[0]->GetYaxis()->SetRangeUser(m2/10 , M2*2);
    //h2[0]->SetStats(0);
    for(int i=0; i<6; i++) h2[i]->Draw("HIST same");
    legend2->Draw();

    if(stamp==1){
        int Cartella= system("mkdir -p FAMU");
        c->SaveAs("FAMU/theta.png");
        c->SaveAs("FAMU/theta.pdf");
        c->SaveAs("FAMU/theta.root");

        c1->SaveAs("FAMU/phi.png");
        c1->SaveAs("FAMU/phi.pdf");
        c1->SaveAs("FAMU/phi.root");
        
        c2->SaveAs("FAMU/p1.png");
        c2->SaveAs("FAMU/p1.pdf");
        c2->SaveAs("FAMU/p1.root");
        
        c3->SaveAs("FAMU/p2.png");
        c3->SaveAs("FAMU/p2.pdf");
        c3->SaveAs("FAMU/p2.root");
        
        c4->SaveAs("FAMU/qp1.png");
        c4->SaveAs("FAMU/qp1.pdf");
        c4->SaveAs("FAMU/qp1.root");
        
        c5->SaveAs("FAMU/qp2.png");
        c5->SaveAs("FAMU/qp2.pdf");
        c5->SaveAs("FAMU/qp2.root");
    }


    return;
}