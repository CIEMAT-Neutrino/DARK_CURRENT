{  // macro to Fit Histograms of SiPM dT distributions in Log dT scale
   // using Exponential distribution transformed into Log(dT)  x-axis
   //+using Negative Binomial derivation transformed to Log(dT) x-axis
   // by P.Filip. 3.Sept. 2020, odd. 33 FZU //

   TFile* TA = new TFile("sipm_DC_filt.root","READ"); // open file
   cout << "*** Attempting to Read SiPM 'hDt' histogram from File: " << TA->GetName() << endl;
   TH1F* dTimeLogAC = (TH1F*)TA->Get("hDt"); // get Histogram

  // TH1D *dTimeLogAC =0;
  // dTimeLogAC = hcas;
  // .x wFit_SiPM.C 
  //
  //
   dTimeLogAC->Draw(); dTimeLogAC->SetLineStyle(1);  // Show simulated distribution

   float sRATE = 0.2;      // guess Fit position of Peak:  Log(<dT>) = Log(1/sRATE)  
   float sMax = 410*1000;  // guess Peak Maximum Value of Histogram being fitted
   float fRATE = 0.0; // fitted Rate Value [Hz]

   float fRATE1 = 0.0; // resulting fitted Rate A, B values [Hz]
   float fRATE2 = 0.0;
   float fMax1  = 0.0; // resulting fitted Peak A, B magnitudes
   float fMax2  = 0.0;

   float eRATE1 = 0.0; // resulting fitted Rates error estimate
   float eRATE2 = 0.0;
   float eMax1  = 0.0; // fitted Peak magnitudes error estimate
   float eMax2  = 0.0;

/// FITTING Peak A
    sRATE = 1.0;    // guessed Rate
    sMax  = 100; // Magnitude
    sRATE = TMath::Power(10.0,-1*(0.5));    // guessed Rate
    sMax  = 100; // Magnitude

   TF1* fG1 = new TF1("fG1","[1]*[0]*log(10)*(10**x)*exp(-[0]*(10**x))");
      fG1->SetParameter(0,sRATE); // init fit value = Rate [Hz]
      fG1->SetParameter(1,sMax);  // init fit value = Height of the Peak Histogram
      fG1->SetRange(-1,2);    // RANGE-A of Fit Function [ Log(dT) units ]

      fG1->SetLineWidth(3); fG1->SetLineColor(3);

   cout << endl << "################  Fiting Peak-A at guessed position Log<dT> = " << flush;
   cout << TMath::Log10(1/sRATE) << " ###########################" << endl; 

   dTimeLogAC->Fit("fG1","R"); // FIT now, using the RANGEa

   fRATE1 = fG1->GetParameter(0); // remember Fit results A
   fMax1  = fG1->GetParameter(1);

   cout << " -------> Fitted Rate A = " << fRATE1 << " Hz,   < dT > = " << 1/fRATE1 << " s." << flush;
   cout << "   Log<dT> = " << TMath::Log10(1/fRATE1) << " <---------" << endl;  
   cout << "################################################################################################" << endl;
   cout << " -------> Fitted DCR =  " << fRATE1/36.0*1000 << " mHz/mm2 " << endl;
   cout << "################################################################################################" << endl;


/// FITTING Peak B using Negative Binomial Distribution - Prior for N=0 probability
   sRATE = 1000.0;  // guessed Rate
//    sMax  = 100; // Magnitude      
   sMax  = 100; // Magnitude      
   float rShape = 120;  // initial overdispersion shape parameter

   TF1* fG2 = new TF1("fG2","[1]*[0]*(([2]-1)/[2])*(10**x)*log(10)*(1+([0]*(10**x))/[2])**(-[2])");
      fG2->SetParameter(0,sRATE);  // init Peak average Rate = 200 Hz;   Log[<dT>] = Log[1/Rate] 
      fG2->SetParameter(1,sMax); // guess Maximum of Peak Histogram
      fG2->SetParameter(2,rShape);  // <----- Shape parameter "r"
      fG2->SetRange(-5.0,-1.0);    // set RANGE-B
//      fG2->SetRange(-4.5,-1.0);    // set RANGE-B

      fG2->SetLineWidth(3); fG2->SetLineColor(5);

   cout << endl << "################  Fiting Peak-B at guessed position Log<dT> = " << flush;
   cout << TMath::Log10(1/sRATE) << " ##########################" << endl;

   dTimeLogAC->Fit("fG2","R"); // FIT using the RANGEb

   fRATE2 = fG2->GetParameter(0);  // remember Fit results B
   fMax2  = fG2->GetParameter(1);
   float r = fG2->GetParameter(2);

   cout << " -------> Fitted Rate B = " << fRATE2 << " Hz,   < dT > = " << 1/fRATE2 << " s." << flush;
   cout << "   Log<dT> = " << TMath::Log10(1/fRATE2) << " <---------" << flush;
   cout << " shape  r = " << r << endl;

/// ********  FITTING Peak A + B simultaneously, Overlap Possible ***************

   TF1* fG22 = new TF1("fG22",  // defining the Fit Function with Two Peaks here
   "[1]*[0]*log(10)*(10**x)*exp(-[0]*(10**x)) + [3]*[4]*(([4]-1)/[4])*(10**x)*log(10)*(1+([2]*(10**x))/[4])**(-[4])");

   fG22->SetParameter(0,fRATE1);  // init fit Rate using partial Fits A, B
   fG22->SetParameter(2,fRATE2);
   fG22->SetParameter(1,fMax1);   // init Peak Magnitudes using Fits A, B
   fG22->SetParameter(3,fMax2);
   fG22->SetParameter(4,1.7); // Guessed Shape Parameter = SHOULD BE 1.7 for real SiPM data

   fG22->SetRange(-4.5, 2.0); // using FULL range ( containing Peak A + Peak B )

   cout << endl;
   cout << "################################################################################################" << endl;
   cout << "################  Fiting Peaks A+B using positions: Log<dT> = " << flush;           
   cout << TMath::Log10(1/fRATE1) << "  and  " << TMath::Log10(1/fRATE2) << " ##########" << endl;

   fG22->SetLineWidth(3); fG22->SetLineColor(2); fG22->SetLineStyle(7);

   dTimeLogAC->Fit("fG22","R"); // FIT now, using the RANGEa
   
   fRATE1 = fG22->GetParameter(0); eRATE1 = fG22->GetParError(0);
   fRATE2 = fG22->GetParameter(2); eRATE2 = fG22->GetParError(2);
   fMax1  = fG22->GetParameter(1); eMax1  = fG22->GetParError(1);
   fMax2  = fG22->GetParameter(3); eMax2  = fG22->GetParError(3);

   cout << " -------> Fitted new Rate A = " << fRATE1 << " Hz,   < dT > = " << 1/fRATE1 << " s." << flush;
   cout << "   Log<dT> = " << TMath::Log10(1/fRATE1) << " <---------" << endl;
   cout << "           Estimated err: +/- " << eRATE1 << endl;
   cout << "                Magnitude A = " << fMax1 << " +/- " << eMax1 << " counts " << endl;
   cout << " -------> Fitted new Rate B = " << fRATE2 << " Hz,   < dT > = " << 1/fRATE2 << " s." << flush;
   cout << "   Log<dT> = " << TMath::Log10(1/fRATE2) << " <---------" << endl;
   cout << "           Estimated err: +/- " << eRATE2 << endl;
   cout << "                Magnitude B = " << fMax2 << " +/- " << eMax2 << " counts." << endl;
   cout << "################################################################################################" << endl;
   cout << "   shape Parameter 'r' = " << fG22->GetParameter(4) << " +/- " << fG22->GetParError(4) << endl;
   cout << "################################################################################################" << endl;
   cout << " -------> Fitted DCR =  " << fRATE1/36.0*1000 << " mHz/mm2 " << endl;
   cout << "################################################################################################" << endl;
   cout << endl;

      fG1->Draw("SAME");
      fG2->Draw("SAME");

      fG22->Draw("SAME"); // the Last Fit is drawn anyway
   
   TLegend *legend = new TLegend(0,1.0,0.3,0.8);
   legend->AddEntry(dTimeLogAC," SiPM Log(dT)","l");
   legend->AddEntry(fG1," Fitted single Peak A","l");
   legend->AddEntry(fG2," Fitted single Peak B","l");
   legend->AddEntry(fG22," Fitted Peaks [ A+B ]","l");
   legend->Draw();

   TText* T1 = new TText( 0.8,420*1000,"A"); T1->Draw();
   TText* T2 = new TText(-1.6,700*1000,"B"); T2->Draw();
   TText* T3 = new TText(-5.9,120*1000,"C"); T3->Draw();

   char buffer [250];
   int n; 
   float a = fRATE1/36.0*1000;
   float b = fRATE1;
   float c = TMath::Log10(1/fRATE1);

   n=sprintf (buffer, "DCR = %5.3f mHz/mm2 (%5.3fHz/log= %5.3f)", a, b, c);
   printf ("[%s] \n",buffer);
   TText* T5 = new TText(-7.5,20,buffer); T5->Draw();

   gPad->Update();
   TPaveStats *st8 = (TPaveStats*)dTimeLogAC->FindObject("stats");
   st8->SetOptStat(111111);
   st8->SetOptFit(1);
   st8->SetY1NDC(0.5);
   gPad->Update();
}
