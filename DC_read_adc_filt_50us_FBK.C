{
#include "Riostream.h"
#include "TObject.h"
#include <vector>

  // using namespace std;

  ///////////////////////////////////////////////////

  TFile *fil2 = new TFile("sipm_DC_filt.root", "RECREATE");

  ////// Booking Histograms //////

  int Nunits = 25002; // Nsamples for 50us
  float Min = 0.;
  double DeltaT = 2.;
  int Nsegs = 50;
  int Nevents = 1;
  // int Nevents = 10.;

  int nsamples = Nunits;

  float Max = Min + Nunits * DeltaT;

  TH1F *hAux;
  TH1F *hW_av[12000]; // Nsegs*Nevents
  TH1F *hW_adc[12000];
  TH1F *hA;
  TH1F *hQ;

  TH1F *hrate = new TH1F("hrate", ";DCR (mHz/mm^{2});Entries", 500, 0., 1000.);

  // A histogram for every segment
  for (int ie = 0; ie < Nevents; ie++)
  {
    for (int i = 0; i < Nsegs; i++)
    {
      hW_av[(Nsegs * ie) + i] = new TH1F(Form("hW_Ev%dSeg%d", ie, i), "", nsamples, Min, Max);
      hW_adc[(Nsegs * ie) + i] = new TH1F(Form("hWadc_Ev%dSeg%d", ie, i), "", nsamples, Min, Max);
    }
  }

  // Logaritmic histogram for Dt vs Amplitude
  const Int_t nbins = 1000;
  const Int_t nbins2 = 100;
  Double_t xmin = 1.e-9;
  Double_t xmax = 10.;
  Double_t logxmin = TMath::Log10(xmin);
  Double_t logxmax = TMath::Log10(xmax);
  Double_t binwidth = (logxmax - logxmin) / nbins;
  Double_t xbins[nbins + 1];
  xbins[0] = xmin;
  for (Int_t i = 1; i <= nbins; i++)
  {
    xbins[i] = xmin + TMath::Power(10, logxmin + i * binwidth);
  }
  TH2F *hDelta = new TH2F("hDelta", ";#DeltaT (s);Amplitude (V)", nbins, xbins, 500., 0., 0.1);

  Double_t binwidth2 = (logxmax - logxmin) / nbins2;
  Double_t xbins2[nbins2 + 1];
  xbins2[0] = xmin;
  for (Int_t i = 1; i <= nbins2; i++)
  {
    xbins2[i] = xmin + TMath::Power(10, logxmin + i * binwidth2);
  }
  TH1F *hDt = new TH1F("hDt", ";#DeltaT (s);", nbins2, xbins2);

  ///////////////////////////////

  TString file;

  TString Line1, Line2, Line3, Line4;
  TString C1, C2;
  Double_t Dt[50]; // Nsegs
  Double_t tim, adc;
  char c;

  double peakA[11000]; // Vectors with pulse information
  double peakT[11000];
  double Delta[11000];
  int p = 0; // number of pulses

  //  float alpha= 25./(25.+2.); // high-pass filter
  float alpha = 150. / (150. + 2.); // high-pass filter

  float rate[200]; // DC rate per event
  float RT = 0.;

  //  for (int i=0;i<Nevents; i++)
  for (int i = 0; i < 1; i++)
  {

    file = "/pc/choozdsk01/palomare/SiPM/TT_FBK_02/DC/OV25-sec-50-50us--";
    // file = "/pc/choozdsk01/palomare/SiPM/DC_SC_FBKTT/OV3--";
    //      if ( i==35 ) continue;

    if (i < 10)
      file += 0000;
    if (i < 100)
      file += 000;
    if (i < 1000)
      file += 00;
    if (i < 10000)
      file += 0;

    file += i;
    file += ".txt";

    cout << "Opening File " << i << "  " << file << endl;
    ifstream fd1;
    fd1.open(file);

    fd1 >> Line1;
    fd1 >> Line2;
    fd1 >> Line3;

    cout << Line1 << " // " << Line2 << " // " << Line3 << endl;

    // Reading segment information

    int Nes = 0;
    for (int ij = 0; ij < Nsegs; ij++)
    {
      fd1 >> C1 >> C2;
      TString C3(C2(9, 17));
      Dt[ij] = C3.Atof();
      cout << C1 << "  " << Dt[ij] << endl;
    }

    //      rate[i] = Nsegs/Dt[Nsegs-1];
    // RT+=rate[i];

    fd1 >> Line4;
    cout << Line4 << endl;

    double sum = 0.;
    int is = 0;
    int l = 0;

    int idc = 0;
    double adc_min = 0.;
    double time_min = 0.;
    double sum_time = 0.;
    double adc_diff = 0.;
    double time_diff = 0.;

    double f = 0.;
    double f_old = 0.;
    double sum_old = 0.;

    for (int ij = 0; ij < Nsegs; ij++)
    {

      if (hAux)
      {
        delete hAux;
        hAux = NULL;
      }
      hAux = new TH1F("hAux", "", nsamples, Min, Max);

      sum_time = 0.;
      Double_t mean_ped_i = 0.;
      Double_t mean_ped_f = 0.;
      int j = 0;

      double Vadc[25002];

      int Nnoise = 0;
      while (j < Nunits)
      {
        fd1 >> tim >> c >> adc;

        hAux->Fill(DeltaT * (j + 0.5), adc);

        //	      if ( adc>0.002 ) Nnoise++;

        Vadc[j] = adc;

        j++;
      }

      int maximumS1 = hAux->GetMaximumBin();
      Double_t MaxS1 = hAux->GetBinContent(maximumS1);
      if (MaxS1 > 0.003)
        cout << "maximumS1 " << MaxS1 << endl;
      if (MaxS1 < 0.003)
        Nes++;
      //	  cout << "Nnoise " << Nnoise << endl;

      j = 0;
      while (j < Nunits && MaxS1 < 0.003)
      //	  while (j<Nunits)
      {

        // Pedestal info
        if (j < 20)
          mean_ped_i += Vadc[j];
        if (j == 20)
        {
          mean_ped_i /= 20.;
        }
        if (j > 480 && j < 500)
          mean_ped_f += Vadc[j];
        if (j == 500)
          mean_ped_f /= 20.;

        // High-pass filtering
        if (f == 0)
        {
          f = Vadc[j];
          f_old = f;
        }
        else
        {
          f = alpha * f_old + alpha * (Vadc[j] - sum_old);
          f_old = f;
          sum_old = Vadc[j];
        }

        hW_av[(Nsegs * i) + ij]->Fill(DeltaT * (j + 0.5), f);
        hW_adc[(Nsegs * i) + ij]->Fill(DeltaT * (j + 0.5), Vadc[j]);

        // Searching for pulses
        if ((f < -0.0035 && sum_time < 200.) || (sum_time > 0 && sum_time < 200.)) // Threshold and timing should be tuned
                                                                                   //	      if ( (f<-0.003 && sum_time<200.) || (sum_time>0 && sum_time<100.) ) // Threshold and timing should be tuned
        {

          if (f < adc_min)
          {
            adc_min = f;
            time_min = Dt[ij] * 1.e+9 + DeltaT * (j + 0.5);
          }

          l++;
          sum_time = DeltaT * l;

          //		  l++;
        }
        else if (sum_time >= 200.)
        {

          peakA[p] = adc_min;
          peakT[p] = time_min;
          if (p > 0)
            Delta[p] = time_min - peakT[p - 1];

          if (Delta[p] > 200.)
            idc++;
          cout << ">> p " << p << " Seg " << ij << " A " << peakA[p] << " T " << peakT[p] << " D " << Delta[p] << endl;
          hDelta->Fill(Delta[p] * 1.e-9, -1. * peakA[p]);
          hDt->Fill(Delta[p] * 1.e-9);

          p++;

          sum_time = 0.;
          adc_min = 0.;
          time_min = 0.;
          l = 0;
        }

        j++;
      } // end loop one segment

    } // end loop segments

    rate[i] = Nes / Dt[Nsegs - 1];

    float dcr = rate[i] * 1000. / 36.;
    hrate->Fill(dcr);
    RT += rate[i];

  } // end loop Events

  for (int i = 0; i < Nevents; i++)
    cout << i << "  " << rate[i] << endl;

  RT /= float(Nevents);
  cout << "Total DC Rate " << RT << endl;

  fil2->Write();
  fil2->Close();
}
