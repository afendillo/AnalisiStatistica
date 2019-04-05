#define Pi TMath::Pi()
#define NFibre 32
#define spessore 0.3 //cm
#define lunghezza 10 //cm
#define FAC 0.04
#define Q_e 1.60217662e-19
#define L 210 //cm

using namespace std;
struct posizione3D{double x , y , theta, phi;int piano, fibra, flag=0;};
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
    if (pos.theta!=0 && pos.theta!=Pi/2)
    {
        NewPos.x=pos.x-1*sin(NewPos.phi)*spessore/tan(NewPos.theta);//Asse Z positivo verso l'alto, i muoni vanno verso il basso
        NewPos.y=pos.y-1*cos(NewPos.phi)*spessore/tan(NewPos.theta);
        NewPos.fibra=FibraPos(NewPos);
    }
    else if(pos.theta==Pi/2) NewPos=pos;

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


void FAMU(){

    Reset();
    double I0=1;//cm-2s-1
    double G=1e6;//Gain
    double QE=0.5;//efficienza di conversione
    double Conv = 180/Pi;

    TRandom3* rand= new TRandom3(1);
    posizione3D posPiano1Up, posPiano1D , posPiano2D, posPiano2Up;

    int n=0;

    TGraph2D* graph = new TGraph2D();
    posPiano1Up = CoordPoint1P(rand); 

    graph->SetPoint(n ,posPiano1Up.x,posPiano1Up.y,0.9 );

    posPiano1Up.piano=1;
    posPiano1D= Proiezione(posPiano1Up);
    if (posPiano1D.flag==0) {
        graph->SetPoint(n+1 ,posPiano1D.x,posPiano1D.y,0.61 );
        posPiano2Up =InversionXY(posPiano1D);
        posPiano2Up.piano=2;
        graph->SetPoint(n+2 ,posPiano2Up.x,posPiano2Up.y,0.59 );

        posPiano2D=Proiezione(posPiano2Up);
        posPiano2D.piano=2;
        if(posPiano2D.flag==0){
            graph->SetPoint(n+3 ,posPiano2D.x,posPiano2D.y,0.3 );
        }
    }
    cout<<"################################################################################################################## \n";
    cout<<"\n PIANO "<<posPiano1Up.piano<<" Up \n";
    cout<<"\n X : "<<posPiano1Up.x<<" , Y : "<<posPiano1Up.y<<" Angolo theta: "<< Conv*posPiano1Up.theta<<"° , Angolo phi: "<<Conv*posPiano1Up.phi<<"°"<<endl;
    cout<<"Fibra numero: "<<posPiano1Up.fibra<<" , distanza da PM: "<<lunghezza-posPiano1Up.y<<"\n";

    cout<<"\n PIANO "<<posPiano1D.piano<<" D \n";
    cout<<"\n X : "<<posPiano1D.x<<" , Y : "<<posPiano1D.y<<" Angolo theta: "<<Conv*posPiano1D.theta<<"° , Angolo phi: "<<Conv*posPiano1D.phi<<"° \n";
    if (posPiano1D.flag==0) cout<<"Fibra numero: "<<posPiano1D.fibra<<" , distanza da PM: "<<lunghezza-posPiano1D.y<<"\n";
    else cout<<"Fuori accettanza \n";
    cout<<"------------------------------------------------ \n";
    if (posPiano1D.flag==0) {

        cout<<"\n PIANO "<<posPiano2Up.piano<<" D \n";
        cout<<"\n X : "<<posPiano2Up.x<<" , Y : "<<posPiano2Up.y<<" Angolo theta: "<< Conv*posPiano2Up.theta<<"° , Angolo phi: "<<Conv*posPiano2Up.phi<<"° \n";
        cout<<"Fibra numero: "<<posPiano2Up.fibra<<" , distanza da PM: "<<lunghezza-posPiano2Up.x<<"\n";

        cout<<"\n PIANO "<<posPiano2D.piano<<" D \n";
        cout<<"\n X : "<<posPiano2D.x<<" , Y : "<<posPiano2D.y<<" Angolo theta: "<< Conv*posPiano2D.theta<<"° , Angolo phi: "<<Conv*posPiano2D.phi<<"° \n";
        if (posPiano2D.flag==0) cout<<"Fibra numero: "<<posPiano2D.fibra<<" , distanza da PM: "<<lunghezza-posPiano2D.x<<"\n";
        else cout<<"Fuori accettanza \n";
    }
    graph->Draw("LINE");
    
    // for(int i=0; i<10;i++){
    //     posPiano1Up = CoordPoint1P(rand);
    //     posPiano1Up.piano=1;
    //     posPiano1D= Proiezione(posPiano1Up);

    //     //double *p=CoordinateFibra12Theta(posPiano1Up);
    //     cout<<"################################################################################################################## \n";
    //     cout<<"\n PIANO "<<posPiano1Up.piano<<" Up \n";
    //     cout<<"\n X : "<<posPiano1Up.x<<" , Y : "<<posPiano1Up.y<<" Angolo theta: "<< Conv*posPiano1Up.theta<<"° , Angolo phi: "<<Conv*posPiano1Up.phi<<"°"<<endl;
    //     cout<<"Fibra numero: "<<posPiano1Up.fibra<<" , distanza da PM: "<<lunghezza-posPiano1Up.y<<"\n";

    //     cout<<"\n PIANO "<<posPiano1D.piano<<" D \n";
    //     cout<<"\n X : "<<posPiano1D.x<<" , Y : "<<posPiano1D.y<<" Angolo theta: "<<Conv*posPiano1D.theta<<"° , Angolo phi: "<<Conv*posPiano1D.phi<<"° \n";
    //     if (posPiano1D.flag==0) cout<<"Fibra numero: "<<posPiano1D.fibra<<" , distanza da PM: "<<lunghezza-posPiano1D.y<<"\n";
    //     else cout<<"Fuori accettanza \n";
    //     cout<<"------------------------------------------------ \n";
    //     if (posPiano1D.flag==0) {

    //         //NOTA: x e y si sono invertite perchè il secondo piano è ruotato, quindi per usare le funzione è necessario scambiarle, vengono poi nuovamente
    //         //invertite per mantenere coerenza con la scelta degli assi---->Funzione inversioneXY (viene usata solo per non dover ridefinire la funzione
    //         //FibraPos)
    //         posPiano2Up =InversionXY(posPiano1D);
    //         posPiano2Up.piano=2;

    //         posPiano2D=Proiezione(posPiano2Up);
    //         posPiano2D.piano=2;

    //         cout<<"\n PIANO "<<posPiano2Up.piano<<" D \n";
    //         cout<<"\n X : "<<posPiano2Up.x<<" , Y : "<<posPiano2Up.y<<" Angolo theta: "<< Conv*posPiano2Up.theta<<"° , Angolo phi: "<<Conv*posPiano2Up.phi<<"° \n";
    //         cout<<"Fibra numero: "<<posPiano2Up.fibra<<" , distanza da PM: "<<lunghezza-posPiano2Up.x<<"\n";

    //         cout<<"\n PIANO "<<posPiano2D.piano<<" D \n";
    //         cout<<"\n X : "<<posPiano2D.x<<" , Y : "<<posPiano2D.y<<" Angolo theta: "<< Conv*posPiano2D.theta<<"° , Angolo phi: "<<Conv*posPiano2D.phi<<"° \n";
    //         if (posPiano2D.flag==0) cout<<"Fibra numero: "<<posPiano2D.fibra<<" , distanza da PM: "<<lunghezza-posPiano2D.x<<"\n";
    //         else cout<<"Fuori accettanza \n";
    //     }
    //     delete [] p;
    // }
     return;
}