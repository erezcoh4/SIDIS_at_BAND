#include "NPiEvent.hh"



void NPiEvent::SetOutputTTrees(TTree* outTree, bool piplus){
    outTree->Branch("eventnumber"          ,&evnum                 );
    outTree->Branch("runnum"               ,&runnum                );
    outTree->Branch("inclusive"            ,&inclusive             );
    outTree->Branch("e_E_PCAL"             ,&e_E_PCAL              );
    outTree->Branch("e_E_ECIN"             ,&e_E_ECIN              );
    outTree->Branch("e_E_ECOUT"            ,&e_E_ECOUT             );
    outTree->Branch("e_PCAL_W"             ,&e_PCAL_W              );
    outTree->Branch("e_PCAL_V"             ,&e_PCAL_V              );
    outTree->Branch("e_PCAL_x"             ,&e_PCAL_x              );
    outTree->Branch("e_PCAL_y"             ,&e_PCAL_y              );
    outTree->Branch("e_PCAL_z"             ,&e_PCAL_z              );
    outTree->Branch("e_PCAL_sector"        ,&e_PCAL_sector         );
    outTree->Branch("e_DC_sector"          ,&e_DC_sector           );
    outTree->Branch("e_DC_Chi2N"           ,&e_DC_Chi2N            );
    outTree->Branch("e_DC_x"               ,&e_DC_x                , "e_DC_x[3]/D"         );
    outTree->Branch("e_DC_y"               ,&e_DC_y                , "e_DC_y[3]/D"         );
    outTree->Branch("e_DC_z"               ,&e_DC_z                , "e_DC_z[3]/D"         );
    outTree->Branch("DC_layers"            ,&DC_layers             , "DC_layers[3]/I"      );
    outTree->Branch("e"                    ,&e                     );
    outTree->Branch("Ve"                   ,&Ve                    );
    outTree->Branch("Beam"                 ,&Beam                  );
    outTree->Branch("beam_helicity"        ,&beam_helicity         );
    outTree->Branch("q"                    ,&q                     );
    outTree->Branch("Ebeam"                ,&Ebeam                 );
    outTree->Branch("xB"                   ,&xB                    );
    outTree->Branch("Q2"                   ,&Q2                    );
    outTree->Branch("omega"                ,&omega                 );
    outTree->Branch("W"                    ,&W                     );
    outTree->Branch("Z"                    ,Zpips                  );
    outTree->Branch("y"                    ,&y                     );

    outTree->Branch("EventPassedCuts"      ,&EventPassedCuts       );
    outTree->Branch("ePastCutsInEvent"     ,&ePastCutsInEvent      );
    outTree->Branch("Npips"                ,&Npips                 );
    outTree->Branch("Npims"                ,&Npims                 );
    outTree->Branch("Nelectrons"           ,&Ne                    );
    outTree->Branch("Ngammas"              ,&Ngammas               );
    outTree->Branch("Nprotons"             ,&Np                    );
    outTree->Branch("Nneutrons"            ,&Nn                    );

    if (piplus) {
        outTree->Branch("pi_chi2PID"           ,&pips_chi2PID          , "pi_chi2PID[20]/D"    );
        outTree->Branch("pi_PCAL_x"            ,&pips_PCAL_x           , "pi_PCAL_x[20]/D"     );
        outTree->Branch("pi_PCAL_y"            ,&pips_PCAL_y           , "pi_PCAL_y[20]/D"     );
        outTree->Branch("pi_PCAL_z"            ,&pips_PCAL_z           , "pi_PCAL_z[20]/D"     );
        outTree->Branch("pi_PCAL_sector"       ,&pips_PCAL_sector      , "pi_PCAL_sector[20]/D");
        outTree->Branch("pi_DC_sector"         ,&pips_DC_sector        , "pi_DC_sector[20]/D"  );
        outTree->Branch("pi_Chi2N"             ,&pips_Chi2N            , "pi_Chi2N[20]/D"      );
        outTree->Branch("pi_DC_x"              ,&pips_DC_x             , "pi_DC_x[20][3]/D"    );
        outTree->Branch("pi_DC_y"              ,&pips_DC_y             , "pi_DC_y[20][3]/D"    );
        outTree->Branch("pi_DC_z"              ,&pips_DC_z             , "pi_DC_z[20][3]/D"    );
        outTree->Branch("pi_E_PCAL"            ,&pips_E_PCAL           , "pi_E_PCAL[20]/D"     );
        outTree->Branch("pi_E_ECIN"            ,&pips_E_ECIN           , "pi_E_ECIN[20]/D"     );
        outTree->Branch("pi_E_ECIN"            ,&pips_E_ECIN           , "pi_E_ECIN[20]/D"     );
        outTree->Branch("pi_E_ECOUT"           ,&pips_E_ECOUT          , "pi_E_ECOUT[20]/D"    );

        outTree->Branch("Vpi"                  ,&VpiplusArray          );
        outTree->Branch("pi"                   ,&piplusArray           );

        outTree->Branch("eepipsPastKinematicalCuts",&eepipsPastKinematicalCuts ,"eepipsPastKinematicalCuts[20]/O"  );
        outTree->Branch("piPastCutsInEvent"    ,&pipsPastCutsInEvent   ,"piPastCutsInEvent/O"  );
        outTree->Branch("eepipsPastCutsInEvent",&eepipsPastCutsInEvent ,"eepipsPastCutsInEvent/O"  );

        outTree->Branch("piplus_Px"                ,&piplus_Px              , "piplus_Px[20]/D"    );
        outTree->Branch("piplus_Py"                ,&piplus_Py              , "piplus_Py[20]/D"    );
        outTree->Branch("piplus_Pz"                ,&piplus_Pz              , "piplus_Pz[20]/D"    );
        outTree->Branch("piplus_E"                 ,&piplus_E               , "piplus_E[20]/D"    );
        outTree->Branch("Vpiplus_X"                ,&Vpiplus_X              , "Vpiplus_X[20]/D"    );
        outTree->Branch("Vpiplus_Y"                ,&Vpiplus_Y              , "Vpiplus_Y[20]/D"    );
        outTree->Branch("Vpiplus_Z"                ,&Vpiplus_Z              , "Vpiplus_Z[20]/D"    );
    }
    else {
        outTree->Branch("pi_chi2PID"           ,&pims_chi2PID          , "pi_chi2PID[20]/D"    );
        outTree->Branch("pi_PCAL_x"            ,&pims_PCAL_x           , "pi_PCAL_x[20]/D"     );
        outTree->Branch("pi_PCAL_y"            ,&pims_PCAL_y           , "pi_PCAL_y[20]/D"     );
        outTree->Branch("pi_PCAL_z"            ,&pims_PCAL_z           , "pi_PCAL_z[20]/D"     );
        outTree->Branch("pi_PCAL_sector"       ,&pims_PCAL_sector      , "pi_PCAL_sector[20]/D");
        outTree->Branch("pi_DC_sector"         ,&pims_DC_sector        , "pi_DC_sector[20]/D"  );
        outTree->Branch("pi_Chi2N"             ,&pims_Chi2N            , "pi_Chi2N[20]/D"      );
        outTree->Branch("pi_DC_x"              ,&pims_DC_x             , "pi_DC_x[20][3]/D"    );
        outTree->Branch("pi_DC_y"              ,&pims_DC_y             , "pi_DC_y[20][3]/D"    );
        outTree->Branch("pi_DC_z"              ,&pims_DC_z             , "pi_DC_z[20][3]/D"    );
        outTree->Branch("pi_E_PCAL"            ,&pims_E_PCAL           , "pi_E_PCAL[20]/D"     );
        outTree->Branch("pi_E_ECIN"            ,&pims_E_ECIN           , "pi_E_ECIN[20]/D"     );
        outTree->Branch("pi_E_ECIN"            ,&pims_E_ECIN           , "pi_E_ECIN[20]/D"     );
        outTree->Branch("pi_E_ECOUT"           ,&pims_E_ECOUT          , "pi_E_ECOUT[20]/D"    );

        outTree->Branch("Vpi"                  ,&VpiplusArray          );
        outTree->Branch("pi"                   ,&piplusArray           );

        outTree->Branch("eepimsPastKinematicalCuts",&eepimsPastKinematicalCuts ,"eepimsPastKinematicalCuts[20]/O"  );
        outTree->Branch("piPastCutsInEvent"    ,&pimsPastCutsInEvent   ,"piPastCutsInEvent/O"  );
        outTree->Branch("eepimsPastCutsInEvent",&eepimsPastCutsInEvent ,"eepimsPastCutsInEvent/O"  );

        outTree->Branch("piminus_Px"                ,&piminus_Px              , "piminus_Px[20]/D"    );
        outTree->Branch("piminus_Py"                ,&piminus_Py              , "piminus_Py[20]/D"    );
        outTree->Branch("piminus_Pz"                ,&piminus_Pz              , "piminus_Pz[20]/D"    );
        outTree->Branch("piminus_E"                 ,&piminus_E               , "piminus_E[20]/D"    );
        outTree->Branch("Vpiminus_X"                ,&Vpiminus_X              , "Vpiminus_X[20]/D"    );
        outTree->Branch("Vpiminus_Y"                ,&Vpiminus_Y              , "Vpiminus_Y[20]/D"    );
        outTree->Branch("Vpiminus_Z"                ,&Vpiminus_Z              , "Vpiminus_Z[20]/D"    );
    }

    void NPiEvent::InitializeVariables(){
        e = TLorentzVector(0,0,0,db->GetParticle( 11   )->Mass());
        mc_e = TLorentzVector(0,0,0,db->GetParticle( 11   )->Mass());
        
        
        xB          = Q2        = omega     = -9999;
        xF          = y         = M_X       = -9999;
        e_E_ECIN    = e_E_ECOUT = e_E_PCAL  = -9999;
        e_PCAL_W    = e_PCAL_V              = -9999;
        e_PCAL_x    = e_PCAL_y  = e_PCAL_z  = -9999;
        e_PCAL_sector                       = -9999;
        e_DC_sector = e_DC_Chi2N            = -9999;
        for (int regionIdx=0; regionIdx<3; regionIdx++) {
            e_DC_x[regionIdx]               = -9999;
            e_DC_y[regionIdx]               = -9999;
            e_DC_z[regionIdx]               = -9999;
        }
        Ve                                  = TVector3();
        ePastCutsInEvent                    = false;

        piplus          .clear();
        piminus         .clear();
        Vpiplus         .clear();
        Vpiminus        .clear();
        piplusArray     ->Clear();
        piminusArray    ->Clear();
        VpiplusArray    ->Clear();
        VpiminusArray   ->Clear();
        mc_piplus       ->Clear();
        mc_piminus      ->Clear();
        pipluses        .clear();
        pipluses        .clear();
        piminuses       .clear();
        electrons       .clear();
        neutrons        .clear();
        protons         .clear();
        gammas          .clear();
        for (int piIdx=0; piIdx<NMAXPIONS; piIdx++) {
            pips_chi2PID[piIdx]                         = -9999;
            pips_DC_sector[piIdx]                       = -9999;
            pips_PCAL_sector[piIdx]                     = -9999;
            pips_PCAL_W[piIdx] = pips_PCAL_V[piIdx]     = -9999;
            pips_PCAL_x[piIdx] = pips_PCAL_y[piIdx]     = -9999;
            pips_PCAL_z[piIdx]                          = -9999;
            pips_E_PCAL[piIdx]                          = -9999;
            pips_E_ECIN[piIdx] = pips_E_ECOUT[piIdx]    = -9999;
            
            pims_chi2PID[piIdx]                         = -9999;
            pims_DC_sector[piIdx]                       = -9999;
            pims_PCAL_sector[piIdx]                     = -9999;
            pims_PCAL_W[piIdx] = pims_PCAL_V[piIdx]     = -9999;
            pims_PCAL_x[piIdx] = pims_PCAL_y[piIdx]     = -9999;
            pims_PCAL_z[piIdx]                          = -9999;
            pims_E_PCAL[piIdx]                          = -9999;
            pims_E_ECIN[piIdx] = pims_E_ECOUT[piIdx]    = -9999;
            for (int regionIdx=0; regionIdx<3; regionIdx++) {
                pips_DC_x[piIdx][regionIdx]= pips_DC_y[piIdx][regionIdx]    = -9999;
                pips_DC_z[piIdx][regionIdx]                                 = -9999;
                pims_DC_x[piIdx][regionIdx]= pims_DC_y[piIdx][regionIdx]    = -9999;
                pims_DC_z[piIdx][regionIdx]                                 = -9999;
            }
            piplus  .push_back( TLorentzVector(0,0,0,db->GetParticle( 211 )->Mass()) );
            Vpiplus .push_back( TVector3() );
            pipsPastSelectionCuts[piIdx]                = false;
            eepipsPastKinematicalCuts[piIdx]            = false;
            
            piminus .push_back( TLorentzVector(0,0,0,db->GetParticle( -211 )->Mass()) );
            Vpiminus.push_back( TVector3() );
            pimsPastSelectionCuts[piIdx]                = false;
            eepimsPastKinematicalCuts[piIdx]            = false;
            
            piplus_Px[piIdx]    = piplus_Py[piIdx]  = piplus_Pz[piIdx]  = piplus_E[piIdx]   = -9999;
            piminus_Px[piIdx]   = piminus_Py[piIdx] = piminus_Pz[piIdx] = piminus_E[piIdx]  = -9999;
            Vpiplus_X[piIdx]    = Vpiplus_Y[piIdx]  = Vpiplus_Z[piIdx]  = -9999;
            Vpiminus_X[piIdx]   = Vpiminus_Y[piIdx] = Vpiminus_Z[piIdx] = -9999;
            
        }
        DC_layer                                        = -9999;
        status                                          = 1; // 0 is good...
        
        pipsPastCutsInEvent                             = false;
        eepipsPastCutsInEvent                           = false;
        pimsPastCutsInEvent                             = false;
        eepimsPastCutsInEvent                           = false;
    }
}

void NPiEvent::ComputeKinematics(){
    // compute event kinematics (from e-only information)
    q       = Beam - e;
    Q2      = -q.Mag2();
    omega   = q.E();
    xB      = Q2/(2. * Mp * q.E());
    W2      = Mp2 - Q2 + 2. * omega * Mp;
    W       = sqrt(W2);
    y       = omega / Ebeam;
}