///////////////////////////////////////////

// Create a DATABASE_DFh to sum using USER_DIAMETER

List<scalar> DATABASE_DIAMETERS(34); 
List<scalar> DATABASE_Vt(34); 
DATABASE_DIAMETERS = {0, 0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.2,  1.4, 1.6,  1.8,  2.0,  2.2,  2.4,  2.6,  2.8,  3.0,  3.2,  3.4, 3.6,  3.8,  4.0,  4.2,  4.4,  4.6,  4.8,  5.0,  5.2,  5.4, 5.6,  5.8,  6.0 };
DATABASE_Vt = {0, 1.17, 1.62, 2.06, 2.47, 2.87, 3.27, 3.67, 4.03, 4.64, 5.17, 5.65, 6.09, 6.49, 6.90, 7.27, 7.57, 7.82, 8.06, 8.26, 8.44, 8.60, 8.72, 8.83, 8.92, 8.98, 9.03, 9.07, 9.09, 9.12, 9.14, 9.16, 9.17, 9.18 };  

List<scalar> USER_Vt(phases.size());
List<scalar> USER_DIAMETERS(phases.size());

forAll (phases, i)
{
    USER_DIAMETERS[i] = phases[i];
}

int Vt_index = 0;

// This step is probably made redundant
forAll (phases, i)
{    
    for (label j = 0; j < 34; j++)
    {    
        if ( std::abs(DATABASE_DIAMETERS[j]-USER_DIAMETERS[i]) < scalar(1E-4) )
        {            
            USER_Vt[Vt_index] = DATABASE_Vt[j];
            Vt_index = Vt_index + 1;
        }
    }
}

List<scalar> Vt(phases.size());
Vt = USER_Vt;

List<List<scalar>> DATABASE_DFh(Rh.size());
List<List<scalar>> USER_DFh(Rh.size());
List<List<scalar>> DFh(Rh.size());

forAll (Rh, i)
{    
    scalar Rh_ = Rh[i];
    List<scalar> DATABASE_DF(DATABASE_DIAMETERS.size());
    scalar DATABASE_DFh_SUM = 0;
   
    if (Rh_ > 0)
    {
        // Create DATABASE_DFh
        for (label j = 0; j < 34; j++)
        {
            if (j==0)
            {
                DATABASE_DF[j] = best( (DATABASE_DIAMETERS[j] + (DATABASE_DIAMETERS[j+1]-DATABASE_DIAMETERS[j])/2.0) , Rh_) - 0.0;
            }
            else if (j == phases.size()-1)
            {
                DATABASE_DF[j] = 1 - best( (DATABASE_DIAMETERS[j] - (DATABASE_DIAMETERS[j]*1E3-DATABASE_DIAMETERS[j-1])/2.0) , Rh_);
            }
            else
            {
                DATABASE_DF[j] = best( (DATABASE_DIAMETERS[j] + (DATABASE_DIAMETERS[j+1]-DATABASE_DIAMETERS[j])/2.0) , Rh_) - best( (DATABASE_DIAMETERS[j] - (DATABASE_DIAMETERS[j]-DATABASE_DIAMETERS[j-1])/2.0) , Rh_);
            }
            DATABASE_DFh_SUM = DATABASE_DFh_SUM + (DATABASE_DF[j]*DATABASE_Vt[j]);
        }
        
        DATABASE_DFh[i] = (DATABASE_DF*DATABASE_Vt)/DATABASE_DFh_SUM;

        // Create USER_DFh from DATABASE_DFh by summing discrete DFh until the average of 2 diameters. eg for 0.5mm -> sum DFh from 0 till D=0.75mm (the 0.7mm element)
        scalar USER_DFh_SUM = 0;
        forAll (phases, m)
        {    
            if (m < phases.size()-1)
            {
                USER_DFh_SUM = 0;
                for (label n = 0; n < 34; n++)
                {    
                    if ( DATABASE_DIAMETERS[n] < (USER_DIAMETERS[m]*1E3-USER_DIAMETERS[m-1])/2.0 );
                    {            
                        USER_DFh_SUM = USER_DFh_SUM + DATABASE_DFh[i][n];
                    }
                }
                if (m == 0)
                {
                    USER_DFh[i][m] = USER_DFh_SUM;
                } else {
                    // Due to sum being from 0 till M, we need to subtract the previous M-1 elements to get the interval
                    USER_DFh[i][m] = USER_DFh_SUM - USER_DFh[i][m-1];
                }
            } else if (m == phases.size()-1) {
                // Since summation of all element is 1 and the second last sum is recorded under USER_DFh_SUM, we can just subtract this to get the final interval
                USER_DFh[i][m] = 1 - USER_DFh_SUM;
            }
            Info << nl << "Rh = " << Rh_ << nl;
            Info << "D = " << USER_DIAMETERS[m] << " | fhd = " << USER_DFh[i][m] << nl;
            
        }
    }

    // Write USER_DFh to DFh because DFh is used in calculateCatchRatio.H and I don't want to rename it
    DFh[i] = USER_DFh[i];
    
}
///////////////////////////////////////////
