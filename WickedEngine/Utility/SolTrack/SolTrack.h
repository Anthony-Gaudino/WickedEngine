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

#ifndef __SOLTRACK_H
#define __SOLTRACK_H

#include "stdio.h"
#include "math.h"

///////////////////////////////////////////////////////////////////////////////
// Constants:
///////////////////////////////////////////////////////////////////////////////
#define PI      3.14159265358979323846    // Pi
#define TWO_PI  6.28318530717958647693    // 2 pi
#define MPI     3.14159265358979323846e6  // One Megapi...
#define R2D    57.2957795130823208768     // Radians to degrees conversion factor
#define R2H     3.81971863420548805845    // Radians to hours conversion factor



///////////////////////////////////////////////////////////////////////////////
// Structs:
///////////////////////////////////////////////////////////////////////////////

/// @brief Date and time to compute the Sun's position for, in UT
struct Time {
	int year;
	int month;
	int day;
	int hour;
	int minute;

	double second;
};

/// @brief Location to compute the Sun's position for
struct Location {
	double longitude;
	double latitude;

	double sinLat;
	double cosLat;

	double pressure;
	double temperature;
};

/// @brief Position of the Sun
struct Position {
	double julianDay;
	double tJD;
	double tJC;
	double tJC2;

	double longitude;
	double distance;

	double obliquity;
	double cosObliquity;
	double nutationLon;

	double rightAscension;
	double declination;
	double hourAngle;
	double agst;

	double altitude;
	double altitudeRefract;
	double azimuthRefract;

	double hourAngleRefract;
	double declinationRefract;
};

/// @brief Rise and set data for the Sun
struct RiseSet {
	double riseTime;
	double transitTime;
	double setTime;

	double riseAzimuth;
	double transitAltitude;
	double setAzimuth;
};



///////////////////////////////////////////////////////////////////////////////
// Function prototypes:
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief   Main function to compute the position of the Sun
 *
 * @param [in]    time - Struct containing date and time to compute the
 *     position for, in UT
 * @param [inout] location - Struct containing the geographic location to
 *     compute the position for
 * @param [out]   position - Struct containing the position of the Sun in
 *     horizontal (and equatorial if desired) coordinates
 *
 * @param [in]    useDegrees - Use degrees for input and output angular
 *     variables, rather than radians
 * @param [in]    useNorthEqualsZero - Use the definition where azimuth=0
 *     denotes north, rather than south
 * @param [in]    computeRefrEquatorial - Compute refraction correction for
 *     equatorial coordinates: 0-no, 1-yes
 * @param [in]    computeDistance - Compute distance to the Sun (in AU): 0-no,
 *     1-yes
 *
 * Example Usage:
 * @code SolTrack(time, location, &position, 1,1,0,0); @endcode
 */
void SolTrack(
	struct Time      time,
	struct Location  location,
	struct Position* position,
	int              useDegrees,
	int              useNorthEqualsZero,
	int              computeRefrEquatorial,
	int              computeDistance
);

/**
 * @brief  Compute the Julian Day from the date and time
 *
 * Gregorian calendar only (>~1582).
 *
 * @param [in]  year   - Year of date
 * @param [in]  month  - Month of date
 * @param [in]  day    - Day of date
 * @param [in]  hour   - Hour of time
 * @param [in]  minute - Minute of time
 * @param [in]  second - Second of time
 *
 * @retval JulianDay    Julian day for the given date and time
 */
double computeJulianDay(
	int    year,
	int    month,
	int    day,
	int    hour,
	int    minute,
	double second
);

/**
 * @brief  Compute the ecliptic longitude of the Sun for a given Julian Day
 *
 * Also computes the obliquity of the ecliptic and nutation.
 *
 * @param [in] computeDistance - Compute distance to the Sun (in AU): 0-no,
 *     1-yes
 * @param [inout] position - Position of the Sun
 */
void computeLongitude(
	int              computeDistance,
	struct Position* position
);

/**
 * @brief  Convert ecliptic coordinates to equatorial coordinates
 *
 * This function assumes that the ecliptic latitude = 0.
 *
 * @param [in] longitude    - Ecliptic longitude of the Sun (rad)
 * @param [in] cosObliquity - Cosine of the obliquity of the ecliptic
 *
 * @param [out] rightAscension - Right ascension of the Sun (rad)
 * @param [out] declination    - Declination of the Sun (rad)
 */
void convertEclipticToEquatorial(
	double  longitude,
	double  cosObliquity,
	double* rightAscension,
	double* declination
);

/**
 * @brief  Convert equatorial to horizontal coordinates
 *
 * Also corrects for parallax and atmospheric refraction.
 *
 * @param [in]    location - Geographic location of the observer (rad)
 * @param [inout] position - Position of the Sun (rad)
 */
void convertEquatorialToHorizontal(
	struct Location  location,
	struct Position* position
);

/**
 * @brief  Convert equatorial coordinates to horizontal coordinates
 *
 * @param [in]  sinLat - Sine of the geographic latitude of the observer
 * @param [in]  cosLat - Cosine of the geographic latitude of the observer
 * @param [in]  longitude - Geographic longitude of the observer (rad)
 * @param [in]  rightAscension - Right ascension of the Sun (rad)
 * @param [in]  declination - Declination of the Sun (rad)
 * @param [in]  agst - Apparent Greenwich sidereal time (Greenwich mean
 *     sidereal time, corrected for the equation of the equinoxes)
 *
 * @param [out] azimuth - Azimuth ("wind direction") of the Sun (rad; 0=South)
 * @param [out] sinAlt - Sine of the altitude of the Sun above the horizon
 */
void eq2horiz(
	double  sinLat,
	double  cosLat,
	double  longitude,
	double  rightAscension,
	double  declination,
	double  agst,
	double* azimuth,
	double* sinAlt
);

/**
 * @brief  Convert (refraction-corrected) horizontal coordinates to equatorial
 *     coordinates
 *
 * @param [in]  sinLat - Sine of the geographic latitude of the observer
 * @param [in]  cosLat - Cosine of the geographic latitude of the observer
 * @param [in]  azimuth - Azimuth ("wind direction") of the Sun (rad; 0=South)
 * @param [in]  altitude - Altitude of the Sun above the horizon (rad)
 * @param [out] hourAngle - Hour angle of the Sun (rad; 0=South)
 * @param [out] declination - Declination of the Sun (rad)
 */
void convertHorizontalToEquatorial(
	double  sinLat,
	double  cosLat,
	double  azimuth,
	double  altitude,
	double* hourAngle,
	double* declination
);

/**
 * @brief  Convert the South=0 convention to North=0 convention for azimuth and
 *     hour angle
 *
 * South=0 is the default in celestial astronomy.
 * This makes the angles compatible with the compass/wind directions.
 *
 * @param [inout] azimuth - Azimuth ("wind direction") of the Sun (rad)
 * @param [inout] hourAngle - Hour angle of the Sun (rad)
 * @param [in]    computeRefrEquatorial - Compute refraction correction for
 *     equatorial coordinates
 */
void setNorthToZero(double* azimuth,
	double* hourAngle,
	int     computeRefrEquatorial
);

/**
 * @brief  Convert final results from radians to degrees
 *
 * Not touching intermediate results.
 *
 * @param [inout] longitude - Ecliptical longitude of the Sun (rad->deg)
 * @param [inout] rightAscension - Right ascension of the Sun (rad->deg)
 * @param [inout] declination - Declination of the Sun (rad->deg)
 *
 * @param [inout] altitude - Altitude of the Sun above the horizon (rad->deg)
 * @param [inout] azimuthRefract - Azimuth ("wind direction") of the Sun
 *     (rad->deg), corrected for refraction
 * @param [inout] altitudeRefract - Altitude of the Sun above the horizon
 *     (rad->deg), corrected for refraction
 * @param [inout] hourAngle - Hour angle of the Sun (rad->deg)
 * @param [inout] declinationRefract - Declination of the Sun (rad->deg),
 *     corrected for refraction
 * @param [in]    computeRefrEquatorial - Compute refraction correction for
 *     equatorial coordinates
 */
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
);

/**
 * @brief  My version of the atan2() function - ~39% faster than built-in
 *     (in terms of number of instructions)
 *
 * @param [in] y - Numerator of the fraction to compute the arctangent for
 * @param [in] x - Denominator of the fraction to compute the arctangent for
 *
 * atan2(y,x) = atan(y/x), where the result will be put in the correct quadrant
 *
 * @see  https://en.wikipedia.org/wiki/Atan2#Definition_and_computation
 */
//double STatan2(double y, double x);

/**
 * @brief  Compute rise, transit and set times for the Sun, as well as their
 *     azimuths/altitude
 *
 * @param [in]  time - Struct containing date and time to compute the position
 *     for, in UT
 * @param [in]  location - Struct containing the geographic location to compute
 *     the position for
 * @param [out] position - Struct containing the position of the Sun
 * @param [out] riseSet - Struct containing the Sun's rise, transit and set
 *     data
 * @param [in]  rsAlt - Altitude to return rise/set data for (radians; 0. is
 *     actual rise/set). rsAlt>pi/2: compute transit only
 * @param [in]  useDegrees - Use degrees for input and output angular
 *     variables, rather than radians
 * @param [in]  useNorthEqualsZero - Use the definition where azimuth=0 denotes
 *     north, rather than south
 *
 * Example Usage:
 * @code SolTrack(time, location, &position, 0.0); @endcode
 *
 * @note
 * - if rsAlt == 0.0, actual rise and set times are computed
 * - if rsAlt != 0.0, the routine calculates when alt = rsAlt is reached
 * - returns times, rise/set azimuth and transit altitude in the struct
 *   position
 *
 * @see
 * - subroutine riset() in riset.f90 from libTheSky (libthesky.sf.net) for more
 *   info
 */
void SolTrack_RiseSet(
	struct Time      time,
	struct Location  location,
	struct Position* position,
	struct RiseSet*  riseSet,
	double           sa0,
	int              useDegrees,
	int              useNorthEqualsZero
);

/**
 * @brief  Fold an angle to take a value between 0 and 2pi
 *
 * @param [in] angle - Angle to fold (radians)
 */
double rev(double angle);

/**
 * @brief  Fold an angle to take a value between -pi and pi
 *
 * @param [in] angle - Angle to fold (radians)
 */
double rev2(double angle);
#endif
