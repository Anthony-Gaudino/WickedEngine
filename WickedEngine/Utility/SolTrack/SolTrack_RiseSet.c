/*
  SolTrack: a simple, free, fast and accurate C routine to compute the position
  of the Sun

  Copyright (c) 2014-2019  Marc van der Sluys, Paul van Kan and Jurgen
  Reintjes, Sustainable Energy research group, HAN University of applied
  sciences, Arnhem, The Netherlands

  This file is part of the SolTrack package, see:
  http://soltrack.sourceforge.net SolTrack is derived from libTheSky
  (http://libthesky.sourceforge.net) under the terms of the GPL v.3

  This is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option) any
  later version.

  This software is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
  details.

  You should have received a copy of the GNU Lesser General Public License
  along with this code.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "SolTrack.h"

void SolTrack_RiseSet(
	struct Time      time,
	struct Location  location,
	struct Position* position,
	struct RiseSet*  riseSet,
	double			 rsAlt,
	int				 useDegrees,
	int				 useNorthEqualsZero
) {
	int evi;
	int iter;
	double tmdy[3]  = { 0.0, 0.0, 0.0 };
	double cosH0;
	double h0       = 0.0;
	double th0;
	double dTmdy;
	double accur;
	double ha       = 0.0;
	double alt      = 0.0;
	double azalt[3] = { 0.0, 0.0, 0.0 };

	// Compute refraction-corrected equatorial coordinates
	// (Hour angle, declination): 0-no, 1-yes
	int computeRefrEquatorial = 1;

	// Compute the distance to the Sun in AU: 0-no, 1-yes
	int computeDistance = 0;

	double rsa = -0.8333 / R2D; // Standard altitude for the Sun in radians
	if (fabs(rsAlt) > 1.e-9)
	{
		rsa = rsAlt; // Use a user-specified altitude
	}

	// If the used uses degrees, convert the geographic location to radians:
	struct Location llocation = location; // Create local variable
	if (useDegrees)
	{
		llocation.longitude /= R2D;
		llocation.latitude  /= R2D;
	}

	// Set date and time to midnight UT for the desired day:
	struct Time rsTime;
	rsTime.year  = time.year;
	rsTime.month = time.month;
	rsTime.day   = time.day;

	rsTime.hour   = 0;
	rsTime.minute = 0;
	rsTime.second = 0.0;

	SolTrack(
		rsTime,
		llocation,
		position,
		0, // useDegrees = 0: NEVER use degrees internally!
		useNorthEqualsZero,
		computeRefrEquatorial,
		computeDistance
	);

	double agst0 = position->agst; // AGST for midnight

	int evMax = 2; // Compute transit, rise and set times by default (0-2)
	cosH0 = (sin(rsa) - sin(llocation.latitude) * sin(position->declination))
	      / (           cos(llocation.latitude) * cos(position->declination));
	if (fabs(cosH0) > 1.0) // Body never rises/sets
	{
		evMax = 0; // Compute transit time and altitude only
	}
	else
	{
		h0 = rev(2.0 * acos(cosH0)) / 2.0;
	}

	// Transit time in radians; lon0 > 0 for E
	tmdy[0] = rev(
		position->rightAscension - llocation.longitude - position->agst
	);

	if (evMax > 0)
	{
		tmdy[1] = rev(tmdy[0] - h0); // Rise time in radians
		tmdy[2] = rev(tmdy[0] + h0); // Set  time in radians
	}

	for (evi = 0; evi <= evMax; evi++) // Transit, rise, set
	{
		iter = 0;
		// Accuracy;  1e-5 rad ~ 0.14s. Don't make this smaller than 1e-16
		accur = 1.0e-5;

		dTmdy = INFINITY;
		while (fabs(dTmdy) > accur)
		{
			// Solar day in sidereal days in 2000
			th0 = agst0 + 1.002737909350795 * tmdy[evi];

			// Radians -> seconds - w.r.t. midnight (h=0,m=0)
			rsTime.second = tmdy[evi] * R2H * 3600.0;

			SolTrack(
				rsTime,
				llocation,
				position,
				0, // useDegrees = 0: NEVER use degrees internally!
				useNorthEqualsZero,
				computeRefrEquatorial,
				computeDistance
			);  

			// Hour angle
			ha = rev2(th0 + llocation.longitude - position->rightAscension);

			// Altitude
			alt = asin(
				sin(llocation.latitude) * sin(position->declination) +
				cos(llocation.latitude) * cos(position->declination) * cos(ha)
			);

			// Correction to transit/rise/set times:
			if (evi == 0) // Transit
			{
				dTmdy = -rev2(ha);
			}
			else // Rise/set
			{
				dTmdy = (alt - rsa) / (
					  cos(position->declination)
					* cos(llocation.latitude)
					* sin(ha)
				);
			}
			tmdy[evi] = tmdy[evi] + dTmdy;

			// Print debug output to stdOut:
			//printf(" %4i %2i %2i  %2i %2i %9.3lf    ",
			//	rsTime.year,
			//	rsTime.month,
			//	rsTime.day,
			//	rsTime.hour,
			//	rsTime.minute,
			//	rsTime.second
			//);
			//printf(" %3i %4i   %9.3lf %9.3lf %9.3lf \n",
			//	evi,
			//	iter,
			//	tmdy[evi]*24,
			//	fabs(dTmdy)*24,
			//	accur*24
			//);

			iter += 1;

			if (iter > 30) // while loop doesn't converge
			{ 
				break;
			}
		}  // while(fabs(dTmdy) > accur)

		if (iter > 30) // Convergence failed
		{
			printf("\n  *** WARNING:  riset():  Riset failed to converge: \
				%i %9.3lf  ***\n", evi, rsAlt);
			tmdy[ evi] = -INFINITY;
			azalt[evi] = -INFINITY;
		}
		else // Result converged, store it
		{
			if (evi == 0)
			{
				// Transit altitude
				azalt[evi] = alt;
			}
			else
			{
				// Rise,set hour angle -> azimuth
				azalt[evi] = atan2(sin(ha), (
					  cos(ha)					 * sin(llocation.latitude)
					- tan(position->declination) * cos(llocation.latitude)
				));
			}
		}

		if (tmdy[evi] < 0.0 && fabs(rsAlt) < 1.e-9)
		{
			tmdy[evi]  = -INFINITY;
			azalt[evi] = -INFINITY;
		}
	}  // for-loop evi

	// Set north to zero radians for azimuth if desired:
	if (useNorthEqualsZero) {
		azalt[1] = rev(azalt[1] + PI); // Add PI and fold between 0 and 2pi
		azalt[2] = rev(azalt[2] + PI); // Add PI and fold between 0 and 2pi
	}

	// Convert resulting angles to degrees if desired:
	if (useDegrees) {
		azalt[0] *= R2D; // Transit altitude
		azalt[1] *= R2D; // Rise azimuth
		azalt[2] *= R2D; // Set azimuth
	}

	// Store results:
	riseSet->transitTime = tmdy[0] * R2H; // Transit time - radians -> hours
	riseSet->riseTime    = tmdy[1] * R2H; // Rise time - radians -> hours
	riseSet->setTime     = tmdy[2] * R2H; // Set time - radians -> hours

	riseSet->transitAltitude = azalt[0]; // Transit altitude
	riseSet->riseAzimuth     = azalt[1]; // Rise azimuth
	riseSet->setAzimuth      = azalt[2]; // Set azimuth
}



double rev(double angle) {
	// Add 1e6*PI and use the modulo function to ensure 0 < angle < 2pi
	return fmod(angle + MPI, TWO_PI);
}



double rev2(double angle) {
	// Add 1e6*PI and use the modulo function to ensure -pi < angle < pi
	return fmod(angle + PI + MPI, TWO_PI) - PI;
}
