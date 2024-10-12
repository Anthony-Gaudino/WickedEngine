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

void SolTrack(
	struct Time      time,
	struct Location  location,
	struct Position* position,
	int              useDegrees,
	int              useNorthEqualsZero,
	int              computeRefrEquatorial,
	int              computeDistance
) {
	// If the used uses degrees, convert the geographic location to radians:
	struct Location llocation = location; // Create local variable
	if (useDegrees) {
		llocation.longitude /= R2D;
		llocation.latitude  /= R2D;
	}

	// Compute these once and reuse:
	//---------------------------------------------------------------------------
	llocation.sinLat = sin(llocation.latitude);
	// Cosine of a latitude is always positive or zero
	llocation.cosLat = sqrt(1.0 - llocation.sinLat * llocation.sinLat);

	// Compute the Julian Day from the date and time:
	//---------------------------------------------------------------------------
	position->julianDay = computeJulianDay(
		time.year,
		time.month,
		time.day,
		time.hour,
		time.minute,
		time.second
	);

	// Derived expressions of time:
	//---------------------------------------------------------------------------
	// Time in Julian days since 2000.0
	position->tJD = position->julianDay - 2451545.0;
	// Time in Julian centuries since 2000.0
	position->tJC  = position->tJD / 36525.0;
	position->tJC2 = position->tJC * position->tJC; // T^2

	// Compute the ecliptic longitude of the Sun and the obliquity of the
	// ecliptic:
	//---------------------------------------------------------------------------
	computeLongitude(computeDistance, position);

	// Convert ecliptic coordinates to geocentric equatorial coordinates:
	//---------------------------------------------------------------------------
	convertEclipticToEquatorial(
		position->longitude,
		position->cosObliquity,
		&position->rightAscension,
		&position->declination
	);

	// Convert equatorial coordinates to horizontal coordinates, correcting for
	// parallax and refraction:
	//---------------------------------------------------------------------------
	convertEquatorialToHorizontal(llocation, position);

	// Convert the corrected horizontal coordinates back to equatorial
	// coordinates:
	//---------------------------------------------------------------------------
	if (computeRefrEquatorial)
	{
		convertHorizontalToEquatorial(
			llocation.sinLat,
			llocation.cosLat,
			position->azimuthRefract,
			position->altitudeRefract,
			&position->hourAngleRefract,
			&position->declinationRefract
		);
	}

	// Use the North=0 convention for azimuth and hour angle (default: South = 0)
	// if desired:
	//---------------------------------------------------------------------------
	if (useNorthEqualsZero)
	{
		setNorthToZero(
			&position->azimuthRefract,
			&position->hourAngleRefract,
			computeRefrEquatorial
		);
	}

	// If the user wants degrees, convert final results from radians to degrees:
	//---------------------------------------------------------------------------
	if (useDegrees)
	{
		convertRadiansToDegrees(
			&position->longitude,
			&position->rightAscension,
			&position->declination,
			&position->altitude,
			&position->azimuthRefract,
			&position->altitudeRefract,
			&position->hourAngleRefract,
			&position->declinationRefract,
			computeRefrEquatorial
		);
	}
}



double computeJulianDay(
	int year,
	int month,
	int day,
	int hour,
	int minute,
	double
	second
) {
	if (month <= 2)
	{
		year  -=  1;
		month += 12;
	}

	int tmp1 = (int)floor(year / 100.0);
	int tmp2 = 2 - tmp1 + (int)floor(tmp1 / 4.0);

	double dDay = day + hour   / 24.0
	                  + minute / 1440.0
	                  + second / 86400.0;
	double JD = floor(365.250 * (year + 4716)) + floor(30.60010 * (month + 1))
		+ dDay + tmp2 - 1524.5;

	return JD;
}



void computeLongitude(int computeDistance, struct Position* position)
{
	// Mean longitude
	double l0 = 4.895063168 + 628.331966786 * position->tJC
							+   5.291838e-6 * position->tJC2;
	// Mean anomaly
	double m = 6.240060141 + 628.301955152 * position->tJC
	                       -   2.682571e-6 * position->tJC2;

	// Sun's equation of the centre
	double c = (3.34161088e-2 - 8.40725e-5 * position->tJC - 2.443e-7 * position->tJC2)
	* sin(m) + (3.489437e-4   - 1.76278e-6 * position->tJC) * sin(2 * m);

	// True longitude
	double odot = l0 + c;

	// Nutation, aberration:
	//---------------------------------------------------------------------------
	// Lon. of Moon's mean ascending node
	double omg = 2.1824390725 - 33.7570464271 * position->tJC
	                          +  3.622256e-5  * position->tJC2;
	double dpsi = -8.338601e-5 * sin(omg); // Nutation in longitude
	double dist = 1.0000010178;            // Mean distance to the Sun in AU                                                                           

	if (computeDistance)
	{
		// Eccentricity of the Earth's orbit
		double ecc = 0.016708634 - 0.000042037  * position->tJC
		                         - 0.0000001267 * position->tJC2;
		double nu = m + c; // True anomaly
		// Geocentric distance of the Sun in AU
		dist = dist * (1.0 - ecc * ecc    )
		            / (1.0 + ecc * cos(nu));
	}

	double aber = -9.93087e-5 / dist; // Aberration

	// Obliquity of the ecliptic and nutation - do this here, since we've
	// already computed many of the ingredients:
	//-------------------------------------------------------------------------
	// Mean obliquity of the ecliptic
	double eps0 = 0.409092804222 - 2.26965525e-4 * position->tJC
	                             - 2.86e-9       * position->tJC2;
	double deps = 4.4615e-5 * cos(omg); // Nutation in obliquity

	// Save position parameters:
	//-------------------------------------------------------------------------
	// Apparent geocentric longitude, referred to the true equinox of date
	position->longitude = odot + aber + dpsi;
	while (position->longitude > TWO_PI)
	{
		position->longitude -= TWO_PI;
	}

	position->distance = dist; // Distance (AU)

	position->obliquity = eps0 + deps; // True obliquity of the ecliptic
	position->cosObliquity = cos(position->obliquity); // Need the cosine later on
	position->nutationLon = dpsi; // Nutation in longitude
}



void convertEclipticToEquatorial(
	double  longitude,
	double  cosObliquity,
	double* rightAscension,
	double* declination
) {
	double sinLon = sin(longitude);
	// Sine of the obliquity of the ecliptic will be positive in the
	// foreseeable future
	double sinObl = sqrt(1.0 - cosObliquity * cosObliquity);

	// 0 <= azimuth < 2pi
	*rightAscension = atan2(cosObliquity * sinLon, cos(longitude));
	*declination = asin(sinObl * sinLon);
}



void convertEquatorialToHorizontal(
	struct Location  location,
	struct Position* position
) {
	// Greenwich mean sidereal time
	double gmst = 4.89496121 + 6.300388098985 * position->tJD
	                         + 6.77e-6        * position->tJC2;
	// Correction for equation of the equinoxes -> apparent Greenwich sidereal
	// time
	position->agst = gmst + position->nutationLon * position->cosObliquity;


	double sinAlt = 0.0;
	// Azimuth does not need to be corrected for parallax or refraction, hence
	// store the result in the 'azimuthRefract' variable directly:
	eq2horiz(
		location.sinLat,
		location.cosLat,
		location.longitude,
		position->rightAscension,
		position->declination,
		position->agst,
		&position->azimuthRefract,
		&sinAlt
	);

	double alt = asin(sinAlt); // Altitude of the Sun above the horizon (rad)
	// Cosine of the altitude is always positive or zero
	double cosAlt = sqrt(1.0 - sinAlt * sinAlt);

	// Correct for parallax:
	alt -= 4.2635e-5 * cosAlt; // Horizontal parallax = 8.794" = 4.2635e-5 rad
	position->altitude = alt;

	// Correct for atmospheric refraction:
	//-------------------------------------------------------------------------
	// Refraction correction in altitude
	double dalt = 2.967e-4 / tan(alt + 3.1376e-3 / (alt + 8.92e-2));
	dalt *= location.pressure / 101.0 * 283.0 / location.temperature;
	alt  += dalt;

	// to do: add pressure/temperature dependence
	position->altitudeRefract = alt;
}



void eq2horiz(
	double  sinLat,
	double  cosLat,
	double  longitude,
	double  rightAscension,
	double  declination,
	double  agst,
	double* azimuth,
	double* sinAlt
) {
	double ha = agst + longitude - rightAscension; // Local Hour Angle

	// Some preparation, saves ~29%:
	double sinHa = sin(ha);
	double cosHa = cos(ha);

	double sinDec = sin(declination);
	// Cosine of a declination is always positive or zero
	double cosDec = sqrt(1.0 - sinDec * sinDec);
	double tanDec = sinDec / cosDec;

	// 0 <= azimuth < 2pi
	*azimuth = atan2(sinHa, cosHa * sinLat - tanDec * cosLat);
	// Sine of the altitude above the horizon
	*sinAlt = sinLat * sinDec + cosLat * cosDec * cosHa;
}



void convertHorizontalToEquatorial(
	double  sinLat,
	double  cosLat,
	double  azimuth,
	double  altitude,
	double* hourAngle,
	double* declination
) {
	// Multiply used variables:
	double cosAz = cos(azimuth);
	double sinAz = sin(azimuth); // For symmetry

	double sinAlt = sin(altitude);
	// Cosine of an altitude is always positive or zero
	double cosAlt = sqrt(1.0 - sinAlt * sinAlt);
	double tanAlt = sinAlt / cosAlt;

	// Local Hour Angle: 0 <= hourAngle < 2pi
	*hourAngle = atan2(sinAz, cosAz * sinLat + tanAlt * cosLat);
	// Declination
	*declination = asin(sinLat * sinAlt - cosLat * cosAlt * cosAz);
}



void setNorthToZero(
	double* azimuth,
	double* hourAngle,
	int     computeRefrEquatorial
) {
	*azimuth = *azimuth + PI; // Add PI to set North=0

	if (*azimuth > TWO_PI) // Ensure 0 <= azimuth < 2pi
	{
		*azimuth -= TWO_PI;
	}

	if (computeRefrEquatorial)
	{
		*hourAngle = *hourAngle + PI; // Add PI to set North=0

		if (*hourAngle > TWO_PI) // Ensure 0 <= hour angle < 2pi
		{
			*hourAngle -= TWO_PI;
		}
	}
}



void convertRadiansToDegrees(
	double* longitude,
	double* rightAscension,
	double* declination,
	double* altitude,
	double* azimuthRefract,
	double* altitudeRefract,
	double* hourAngle,
	double* declinationRefract,
	int     computeRefrEquatorial
) {
	*longitude              *= R2D;
	*rightAscension         *= R2D;
	*declination            *= R2D;
					        
	*altitude               *= R2D;
	*azimuthRefract         *= R2D;
	*altitudeRefract        *= R2D;

	if (computeRefrEquatorial)
	{
		*hourAngle          *= R2D;
		*declinationRefract *= R2D;
	}
}



//double STatan2(double y, double x)
//{
//	if (x > 0.0)
//	{
//		return atan(y / x);
//	}
//	else if (x < 0.0)
//	{
//		if (y >= 0.0)
//		{
//			return atan(y / x) + PI;
//		}
//		else // y < 0
//		{
//			return atan(y / x) - PI;
//		}
//	}
//	else // x == 0
//	{
//		if (y > 0.0)
//		{
//			return PI / 2.0;
//		}
//		else if (y < 0.0)
//		{
//			return -PI / 2.0;
//		}
//		else // y == 0
//		{
//			return 0.0;
//		}
//	}
//}
