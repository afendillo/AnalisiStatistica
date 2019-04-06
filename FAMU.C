#define Pi TMath::Pi()
#define NFibre 32
#define spessore 0.3 //cm
#define lunghezza 10 //cm
#define FAC 0.04
#define Q_e 1.60217662e-19
#define L 210 //cm
#define EGamma 197.326e-6*2*Pi/420 //mev nm-->h tagliato
#define I0 1//cm-2s-1
#define Conv 180/Pi
static double G=1e6;//Gain
static double QE=0.5;//efficienza di conversione


using namespace std;
struct posizione3D{double x , y , theta, phi;int piano, fibra, flag=0;};
struct segnale {vector <double> Nf, en, q, flag;};
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
    if (pos.theta!=0)
    {
        if((pos.x/spessore)-(int)(pos.x/spessore)==0 && pos.phi>Pi )
        {
            fibra=ceil(pos.x/spessore)-1;
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
    //origine estremo del piano quadrato-->prima fibra = pos 0
    posizione3D pos;
    pos.x=rdm->Rndm()*(NFibre)*spessore;
    pos.y=rdm->Rndm()*lunghezza;
    double *posTP=Angle(rdm);
    pos.theta=posTP[0];
    pos.phi=posTP[1];
    pos.fibra=FibraPos(pos);
    
    return pos;
}

posizione3D Proiezione(posizione3D pos)
{
    posizione3D NewPos;
    NewPos=pos;
    if (pos.theta!=Pi/2)
    {
        NewPos.x=pos.x-1*sin(NewPos.phi)*spessore*tan(NewPos.theta);//Asse Z positivo verso l'alto, i muoni vanno verso il basso
        NewPos.y=pos.y-1*cos(NewPos.phi)*spessore*tan(NewPos.theta);
        NewPos.fibra=FibraPos(NewPos);
        //scelta coord x=R*sin(phi)*sin(theta), y=R*cos(phi)*sin(theta), z=R*cos(theta)--->le coord vanno messe con segno
    }

    if(NewPos.x>9.6 || NewPos.x<0 || NewPos.y<0 || NewPos.y>10) NewPos.flag=1;
    return NewPos;
}
// double* CoordinateFibra12Theta(posizione3D pos)
// {
//     int NF;
//     double *Coord= new double[4];
//     //Hit per il primo piano
//     if (pos.theta!=0 && pos.theta!=Pi && pos.theta!=Pi/2)
//     {
//         Coord[0]=round(abs(pos.x-0.15)/0.3)+1;//fibra colpita
//     }
//     Coord[1]=lunghezza-pos.y;//posizione y della fibra colpita
//     Coord[2]=pos.theta;
//     Coord[4]=pos.phi;
//     return Coord;
// }


posizione3D InversionXY(posizione3D pos){
    posizione3D NewPos;
    NewPos=pos;
    NewPos.x=pos.y;//inverto per poter calcolare correttamente la posizione della fibra
    NewPos.fibra=FibraPos(NewPos);
    NewPos.x=pos.x;//Li riporto a normalità per poter sfruttare lo stesso calcolo nella proiezione
    return NewPos;
}

double QRaccolta(double x, double E){

    double nGamma=E/EGamma;
    return nGamma*FAC*exp(-x/L)*Q_e*G*QE;
}


TGraph2D* deposito(posizione3D pos1 , posizione3D pos2){

    segnale sig;
    TRandom3* rand = new TRandom3(time(0));
    TGraph2D* graph = new TGraph2D();
    double phi= pos1.phi, yout, dist , y=pos1.y , x=pos1.x, z=0, NewX;
    int NumFibre= abs(pos1.fibra-pos2.fibra);
    if(pos1.piano*pos2.piano==1){
        if(pos2.x>pos1.x){
            graph->SetPoint(0, x , y , z);
            cout<<"X: "<<x<<" Y: " <<y <<" Z: "<<z <<endl;
            for(int i=pos1.fibra; i<=pos2.fibra && x<spessore*NFibre && y<lunghezza; i++){
                NewX=i*spessore*(i<pos2.fibra)+pos2.x*(i==pos2.fibra);
                yout=y+(NewX-x)/tan(phi);
                dist=sqrt(pow(x-NewX, 2)+pow(y-NewX/tan(phi),2));
                (sig.en).push_back(dist*1.8);
                (sig.q).push_back(QRaccolta((rand->Rndm())*(yout-y)+y,dist*1.89));//1.89MeV/cm densità energia media depositata
                (sig.Nf).push_back(i);
                (sig.flag).push_back(pos1.piano);
                y=yout;
                z=z+(NewX-x)/(sin(phi)*tan(pos1.theta));
                x=NewX;
                cout<<"X: "<<x<<" Y: " <<y <<" Z: "<<z <<endl;
                graph->SetPoint(i-pos1.fibra+1, x , y , z);
                
            }

        }
    }
    return graph;
}



void FAMU(){

    Reset();

    TRandom3* rand= new TRandom3(1);
    posizione3D posPiano1Up, posPiano1D , posPiano2D, posPiano2Up;

    // int n=0;

    // TGraph2D* graph = new TGraph2D();
    // posPiano1Up = CoordPoint1P(rand); 

    // graph->SetPoint(n ,posPiano1Up.x,posPiano1Up.y,0.9 );

    // posPiano1Up.piano=1;
    // posPiano1D= Proiezione(posPiano1Up);
    // if (posPiano1D.flag==0) {
    //     graph->SetPoint(n+1 ,posPiano1D.x,posPiano1D.y,0.61 );
    //     posPiano2Up =InversionXY(posPiano1D);
    //     posPiano2Up.piano=2;
    //     graph->SetPoint(n+2 ,posPiano2Up.x,posPiano2Up.y,0.59 );

    //     posPiano2D=Proiezione(posPiano2Up);
    //     posPiano2D.piano=2;
    //     if(posPiano2D.flag==0){
    //         graph->SetPoint(n+3 ,posPiano2D.x,posPiano2D.y,0.3 );
    //     }
    // }
    // cout<<"################################################################################################################## \n";
    // cout<<"\n PIANO "<<posPiano1Up.piano<<" Up \n";
    // cout<<"\n X : "<<posPiano1Up.x<<" , Y : "<<posPiano1Up.y<<" Angolo theta: "<< Conv*posPiano1Up.theta<<"° , Angolo phi: "<<Conv*posPiano1Up.phi<<"°"<<endl;
    // cout<<"Fibra numero: "<<posPiano1Up.fibra<<" , distanza da PM: "<<lunghezza-posPiano1Up.y<<"\n";

    // cout<<"\n PIANO "<<posPiano1D.piano<<" D \n";
    // cout<<"\n X : "<<posPiano1D.x<<" , Y : "<<posPiano1D.y<<" Angolo theta: "<<Conv*posPiano1D.theta<<"° , Angolo phi: "<<Conv*posPiano1D.phi<<"° \n";
    // if (posPiano1D.flag==0) cout<<"Fibra numero: "<<posPiano1D.fibra<<" , distanza da PM: "<<lunghezza-posPiano1D.y<<"\n";
    // else cout<<"Fuori accettanza \n";
    // cout<<"------------------------------------------------ \n";
    // if (posPiano1D.flag==0) {

    //     cout<<"\n PIANO "<<posPiano2Up.piano<<" D \n";
    //     cout<<"\n X : "<<posPiano2Up.x<<" , Y : "<<posPiano2Up.y<<" Angolo theta: "<< Conv*posPiano2Up.theta<<"° , Angolo phi: "<<Conv*posPiano2Up.phi<<"° \n";
    //     cout<<"Fibra numero: "<<posPiano2Up.fibra<<" , distanza da PM: "<<lunghezza-posPiano2Up.x<<"\n";

    //     cout<<"\n PIANO "<<posPiano2D.piano<<" D \n";
    //     cout<<"\n X : "<<posPiano2D.x<<" , Y : "<<posPiano2D.y<<" Angolo theta: "<< Conv*posPiano2D.theta<<"° , Angolo phi: "<<Conv*posPiano2D.phi<<"° \n";
    //     if (posPiano2D.flag==0) cout<<"Fibra numero: "<<posPiano2D.fibra<<" , distanza da PM: "<<lunghezza-posPiano2D.x<<"\n";
    //     else cout<<"Fuori accettanza \n";
    // }
    // graph->Draw("LINE");
    
    for(int i=0; i<1;i++){
        posPiano1Up = CoordPoint1P(rand);
        posPiano1Up.piano=1;
        posPiano1D= Proiezione(posPiano1Up);

        auto graph = deposito(posPiano1Up, posPiano1D);

        //cout<<(posPiano1Up.phi>0 && posPiano1Up.phi<Pi)<<" E "<<(posPiano1D.fibra>posPiano1Up.fibra)<<" E "<<(posPiano1D.x>posPiano1Up.x)<<endl;
        //double *p=CoordinateFibra12Theta(posPiano1Up);
        cout<<"################################################################################################################## \n";
        cout<<"\n PIANO "<<posPiano1Up.piano<<" Up \n";
        cout<<"\n X : "<<posPiano1Up.x<<" , Y : "<<posPiano1Up.y<<" Angolo theta: "<< Conv*posPiano1Up.theta<<"° , Angolo phi: "<<Conv*posPiano1Up.phi<<"°"<<endl;
        cout<<"Fibra numero: "<<posPiano1Up.fibra<<" , distanza da PM: "<<lunghezza-posPiano1Up.y<<"\n";

        cout<<"\n PIANO "<<posPiano1D.piano<<" D \n";
        cout<<"\n X : "<<posPiano1D.x<<" , Y : "<<posPiano1D.y<<" Angolo theta: "<<Conv*posPiano1D.theta<<"° , Angolo phi: "<<Conv*posPiano1D.phi<<"° \n";
        if (posPiano1D.flag==0) cout<<"Fibra numero: "<<posPiano1D.fibra<<" , distanza da PM: "<<lunghezza-posPiano1D.y<<"\n";
        else cout<<"Fuori accettanza \n";
        cout<<"------------------------------------------------ \n";
        // double R=sqrt(pow(posPiano1D.x-posPiano1Up.x,2)+pow(posPiano1D.y-posPiano1Up.y,2)+0.3*0.3);
        // cout<<"R: "<<R<<endl;
        // cout<<"X: "<<sqrt(pow(posPiano1Up.x,2)+pow(posPiano1Up.y,2))*sin(posPiano1Up.phi)*sin(posPiano1Up.theta)<<endl;
        // cout<<"y: "<<sqrt(pow(posPiano1Up.x,2)+pow(posPiano1Up.y,2))*sin(posPiano1Up.theta)<<endl;
        // cout<<"z: "<<R*cos(posPiano1Up.theta)<<endl;
        graph->SetTitle("Spost; X; Y; Z ");
        graph->SetMarkerStyle(20);
        graph->Draw("LINE PL");
        // if (posPiano1D.flag==0) {

        //     //NOTA: x e y si sono invertite perchè il secondo piano è ruotato, quindi per usare le funzione è necessario scambiarle, vengono poi nuovamente
        //     //invertite per mantenere coerenza con la scelta degli assi---->Funzione inversioneXY (viene usata solo per non dover ridefinire la funzione
        //     //FibraPos)
        //     posPiano2Up =InversionXY(posPiano1D);
        //     posPiano2Up.piano=2;

        //     posPiano2D=Proiezione(posPiano2Up);
        //     posPiano2D.piano=2;

        //     cout<<"\n PIANO "<<posPiano2Up.piano<<" D \n";
        //     cout<<"\n X : "<<posPiano2Up.x<<" , Y : "<<posPiano2Up.y<<" Angolo theta: "<< Conv*posPiano2Up.theta<<"° , Angolo phi: "<<Conv*posPiano2Up.phi<<"° \n";
        //     cout<<"Fibra numero: "<<posPiano2Up.fibra<<" , distanza da PM: "<<lunghezza-posPiano2Up.x<<"\n";

        //     cout<<"\n PIANO "<<posPiano2D.piano<<" D \n";
        //     cout<<"\n X : "<<posPiano2D.x<<" , Y : "<<posPiano2D.y<<" Angolo theta: "<< Conv*posPiano2D.theta<<"° , Angolo phi: "<<Conv*posPiano2D.phi<<"° \n";
        //     if (posPiano2D.flag==0) cout<<"Fibra numero: "<<posPiano2D.fibra<<" , distanza da PM: "<<lunghezza-posPiano2D.x<<"\n";
        //     else cout<<"Fuori accettanza \n";
        // }
        //delete [] p;
    }
     return;
}