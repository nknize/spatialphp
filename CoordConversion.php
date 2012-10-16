<?php
    
/**
 * Error handling logic
 */
function assertionHandler($file, $line, $code, $desc = null) {
    print("Assertion failed at " . $file . " " . $line . " " . $code);
    if($desc)
        print(": " . $desc);
    print("<br>");
}

assert_options(ASSERT_CALLBACK, 'assertionHandler');

class LatZones {
    private $letters = array('A', 'C', 'D', 'E', 'F', 'G', 'H', 'J', 'K',
             'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Z');
    private $degrees = array(-90, -84, -72, -64, -56, -48, -40, -32, -24, -16,
             -8, 0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 84);
    private $negLetters = array('A', 'C', 'D', 'E', 'F', 'G', 'H', 'J', 'K', 'L', 'M');
    private $negDegrees = array(-90, -84, -72, -64, -56, -48, -40, -32, -24, -16, -8);
    private $posLetters = array('N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Z');
    private $posDegrees = array(0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 84);
    private $arrayLength = 22;

    public function getLatZoneDegree($letter) {
        assert('!empty($letter) && isset($letter)', 'input variable is empty');
        $ltr = $letter{0};
        for($i=0; $i<$this->arrayLength; $i++) {
            if($this->letters[$i] == $ltr)
                return $this->degrees[$i];
        }
        return -100;
    }

    public function getLatZone($latitude) {
        assert('!empty($latitude) && isset($latitude)', 'input variable is empty');
        $latIndex = -2;
        $lat = intval($latitude);
        if($lat >= 0) {
           $len = count($this->posLetters);
           for($i=0; $i<$len; $i++) {
               if($lat == $this->posDegrees[$i]) {
                   $latIndex = $i;
                   break;
               }
               if($lat > $this->posDegrees[$i]) {
                   continue;
               } else {
                   $latIndex = $i-1;
                   break;
               }
           }
        } else {
            $len = count($this->negLetters);
            for($i=0; $i<$len; $i++) {
                if($lat == $this->negDegrees[$i]) {
                    $latIndex = $i;
                    break;
                }
                if($lat<$this->negDegrees[$i]) {
                    $latIndex = $i-1;
                    break;
                } else {
                    continue;
                }
            }
        }
        if($latIndex == -1)
            $latIndex = 0;
        if($lat >= 0) {
            if($latIndex == -2)
                $latIndex = count($this->posLetters) - 1;
            return $this->posLetters[$latIndex] . "";
        } else {
            if($latIndex == -2)
                $latIndex = count($this->negLetters) - 1;
            return $this->negLetters[$latIndex] . "";
        }
    }
}

/**
 * Utility class to convert UTM to Lat/Lon
 */ 
class UTM2LatLon {
    var $easting;
    var $northing;
    var $zone;
    private $southernHem = "ACDEFGHJKLM";

    private $arc;
    private $mu;
    private $ei;
    private $ca;
    private $cb;
    private $cc;
    private $cd;
    private $n0;
    private $r0;
    private $_a1;
    private $dd0;
    private $t0;
    private $Q0;
    private $lof1;
    private $lof2;
    private $lof3;
    private $_a2;
    private $phi1;
    private $fact1;
    private $fact2;
    private $fact3;
    private $fact4;
    private $zoneCM;
    private $_a3;
    private $b = 6356752.314;
    private $a = 6378137;
    private $e = 0.081819191;
    private $e1sq = 0.006739497;
    private $k0 = 0.9996;

    /**
     * Returns the hemisphere of the zone
     */
    protected function getHemisphere($latZone) {
        assert('!empty($latZone) && isset($latZone)', 'input variable is empty');
        if(strrpos($sourthernHem, $latZone))
            return "S";
        else return "N";
    }

    /**
     * Converts input utm to lat lon value
     */
    public function convertUTMToLatLon($utm) {
        assert('!empty($utm) && isset($utm)', 'input variable is empty');
        $latLon = array(0.0, 0.0);
        $utmArr = split(" ", $utm);
        $this->zone = intval($utmArr[0]);
        $latZone = intval($utmArr[1]);
        $this->easting = doubleval($utmArr[2]);
        $this->northing = doubleval($utmArr[3]);
        $hemisphere = $this->getHemisphere($latZone);

        if($hemisphere=="S")
            $this->northing = 10000000 - $this->northing;
        $this->setVariables();
        $latitude = 180 * ($this->phi1 - $this->fact1 * ($this->fact2 + $this->fact3 + $this->fact4)) / pi();
        if($this->zone > 0)
            $this->zoneCM = 6.0 * $this->zone - 183.0;
        else
            $this->zoneCM = 3.0;
        $longitude = $this->zoneCM - $this->_a3;
        if($hemisphere=="S")
            $latitude = -$latitude;
        return array($latitude, $longitude);
    }

    /**
     * Set the independent / dependent variables
     */
    protected function setVariables() {
        $this->arc = $this->northing / $this->k0;
        $this->mu = $this->arc / 
                ($this->a * (1-pow($this->e, 2.0)/4.0-3.0*pow($this->e, 4.0) / 64.0 - 5.0 * pow($this->e, 6.0) / 256.0));
        $this->ei = (1.0 - pow((1.0 - $this->e * $this->e), (1.0 / 2.0))) /
                    (1.0 + pow((1.0 - $this->e * $this->e), (1.0 / 2.0)));
        $this->ca = 3.0 * $this->ei / 2.0 - 27.0 * pow($this->ei, 3.0) / 32.0;
        $this->cb = 21.0 * pow($this->ei, 2.0) / 16.0 - 55.0 * pow($this->ei, 4.0) / 32.0;
        $this->cc = 151.0 * pow($this->ei, 3.0) / 96.0;
        $this->cd = 1097.0 * pow($this->ei, 4.0) / 512.0;
        $this->phi1 = $this->mu + $this->ca * sin(2.0 * $this->mu) + $this->cb * sin(4.0*$this->mu) + 
                      $this->cc * sin(6.0 * $this->mu) + $this->cd * sin(8.0 * $this->mu);
        $this->n0 = $this->a / pow((1.0 - pow(($this->e * sin($this->phi1)), 2.0)), (1.0/2.0));
        $this->r0 = $this->a * (1.0 - $this->e * $this->e) / pow((1.0 - pow(($this->e * sin($this->phi1)), 2.0)), (3.0/2.0));
        $this->fact1 = $this->n0 * tan($this->phi1)/$this->r0;
        $this->_a1 = 500000 - $this->easting;
        $this->dd0 = $this->_a1 / ($this->n0 * $this->k0);
        $this->fact2 = $this->dd0 * $this->dd0 / 2.0;
        $this->t0 = pow(tan($this->phi1), 2.0);
        $this->Q0 = $thise1sq * pow(cos($this->phi1), 2.0);
        $this->fact3 = (5.0+3.0*$this->t0+1.0*$this->Q0-4.0*$this->Q0*$this->Q0-9.0*$this->e1sq)*pow($this->dd0,4.0)/24.0;
        $this->fact4 = (61.0+90.0*$this->t0+298.0*$this->Q0+45.0*$this->t0*$this->t0-
                        252.0*$this->e1sq-3.0*$this->Q0*$this->Q0)*pow($this->dd0,6.0)/720.0;
        $this->lof1 = $this->_a1 / ($this->n0*$this->k0);
        $this->lof2 = (1.0 + 2.0 * $this->t0 + $this->Q0) * pow($this->dd0, 3.0) / 6.0;
        $this->lof3 = (5.0 - 2.0 * $this->Q0 + 28.0 * $this->t0 - 3.0 * pow($this->Q0, 2.0) + 8.0 * $this->e1sq + 
                      24.0 * pow($this->t0, 2.0));
        $this->_a2 = ($this->lof1 - $this->lof2 + $this->lof3) / cos($this->phi1);
        $this->_a3 = $this->_a2 * 180.0 / pi();
    }
}

class LatLon2UTM {

    // equatorial radius
    private static $equatorialRadius = 6378137;
    // polar radius
    private static $polarRadius = 6356752.314;
    // flattening
    private static $flattening = 0.00335281066474748;// (equatorialRadius-polarRadius)/equatorialRadius;
    // inverse flattening 1/flattening
    private static $inverseFlattening = 298.257223563;// 1/flattening;
    // Mean radius
    private $rm; // = pow($equatorialRadius * $polarRadius, 1 / 2.0);
    // scale factor
    private static $k0 = 0.9996;
    // eccentricity
    private $e; // = sqrt(1.0 - pow($polarRadius / $equatorialRadius, 2.0));
    private $e1sq; // = $e * $e / (1 - $e * $e);
    private $n;// = ($equatorialRadius - $polarRadius) / ($equatorialRadius + $polarRadius);
    // r curv 1
    private $rho = 6368573.744;
    // r curv 2
    private $nu = 6389236.914;
    // Calculate Meridional Arc Length
    // Meridional Arc
    private $S = 5103266.421;
    private $A0 = 6367449.146;
    private $B0 = 16038.42955;
    private $C0 = 16.83261333;
    private $D0 = 0.021984404;
    private $E0 = 0.000312705;
    // Calculation Constants
    // Delta Long
    private $p = -0.483084;
    private $sin1; // = ((float)"4.84814e-06");
    // Coefficients for UTM Coordinates
    private $K1 = 5101225.115;
    private $K2 = 3750.291596;
    private $K3 = 1.397608151;
    private $K4 = 214839.3105;
    private $K5 = -2.995382942;
    private $A6; // = ((float)"-1.00541e-07");

    function __construct() {
        // Mean radius
        $this->rm = pow($equatorialRadius * $polarRadius, 1 / 2.0);
        // eccentricity
        $this->e = sqrt(1.0 - pow($polarRadius / $equatorialRadius, 2.0));
        $this->e1sq = $e * $e / (1 - $e * $e);
        $this->n = ($equatorialRadius - $polarRadius) / ($equatorialRadius + $polarRadius);
        $this->sin1 = ((float)"4.84814e-06");
        $this->A6 = ((float)"-1.00541e-07");
    }

    public function validate($latitude, $longitude) {
        assert('!empty($latitude) && isset($latitude)', 'input variable latitude is empty');
        assert('!empty($longitude) && isset($longitude)', 'input variable longitude is empty');
        assert('$latitude >= -90.0 || $latitude <= 90.0 || $longitude >= -180.0 || $longitude < 180.0)', 
               'Legal ranges: latitude [-90,90], longitude [-180, 180)');
    }


    public function convertLatLonToUTM($lat, $lon) {
        assert('!empty($lat) && isset($lat)', 'input variable lat is empty');
        assert('!empty($lon) && isset($lon)', 'input variable lon is empty');

        validate($lat, $lon);
        $this->setVariables($lat, $lon);

        $longZone = getLongZone($lon);
        $latZones = new LatZones();
        $latZone = $latZones.getLatZone($lat);

        $_easting = $this->getEasting();
        $_northing = $this->getNorthing();

        $UTM = $longZone . " " . $latZone . " " . (string)(intval($_easting)) . " " . (string)(intval($_northing));

        return $UTM;
    }

    protected function setVariables($latitude, $longitude) {
        assert('!empty($latitude) && isset($latitutde)', 'input variable lat is empty');
        assert('!empty($longitude) && isset($longitude)', 'input variable lon is empty');
       
        $latitude = degreeToRadian($latitude);
        $this->rho = $this->equitorialRadius / pow(1.0 - pow($this->e * sin($latitude), 2.0), (3.0 / 2.0)); 
        $this->nu = $this->equitorialRadius / pow(1.0 - pow($this->e * sin($latitude), 2.0), (1.0/2.0));

        if($longitude < 0.0)
            $var1 = ((int)((180.0+$longitude)/6.0)) + 1.0;
        else
            $var1 = ((int)($longitude/6.0)) + 31.0;
        $var2 = (6.0*$var1) - 183.0;
        $var3 = $longitude - $var2;
        $this->p = $var3 * 3600.0 / 10000.0;
        $this->S = $this->A0 * $latitude - $this->B0 * sin(2.0 * $latitude) + $this->C0 * sin(4.0 * $latitude) -
                   $this->D0 * sin(6.0 * $latitude) + $this->E0 * sin(8.0 * $latitude);
        $this->K1 = $this->S * $this->k0;
        $this->K2 = $this->nu * sin($latitude) * cos($latitude) * pow($this->sin1, 2.0) * $this->k0 * (100000000.0) / 2.0;
        $this->K3 = ((pow($this->sin1, 4.0) * $this->nu * sin($latitude) * pow(cos($latitude), 3.0)) / 24.0) *
                    (5.0 - pow(tan($latitude), 2.0) + 9.0 * $this->e1sq * pow(cos($latitude), 2.0) + 4.0 * 
                    pow($this->e1sq, 2.0) * pow(cos($latitude), 4.0)) *
                    $this->k0 * (10000000000000000);
        $this->K4 = $this->nu * cos($latitude) * $this->sin1 * $this->k0 * 10000;
        $this->K5 =  pow($this->sin1 * cos($latitude), 3.0) * ($this->nu / 6.0) * 
                     (1.0 - pow(tan($latitude), 2.0) + $this->e1sq * pow(cos($latitude), 2.0)) * $this->k0 * 1000000000000;
        $this->A6 = (pow($this->p * $this->sin1, 6.0) * $this->nu * sin($latitude) * pow(cos($latitude), 5.0) / 720.0) *
                    (61.0 - 58.0 * pow(tan($latitude), 2.0) + pow(tan($latitude), 4.0) + 270.0 *
                    $this->e1sq * pow(cos($latitude), 2.0) - 330.0 * $this->e1sq * pow(sin($latitude), 2.0)) * 
                    $this->k0 * ((float)("1e+24"));
    }

    protected function getLongZone($longitude) {
        assert('!empty($longitude) && isset($longitude)', 'input variable longitude is empty');
        $longZone = 0;
        if($longitude<0.0)
            $longZone = ((180.0 + $longitude) / 6.0) + 1.0;
        else
            $longZone = ($longitude/6.0) + 31.0;
        $val =(string)(intval($longZone));
        if(strlen($val) == 1)
            $val = "0" . $val;
        return $val;
    }

    protected function getNorthing($latitude) {
        assert('!empty($latitude) && isset($latitude)', 'input variable latitude is empty');
        $northing = $this->K1 + $this->K2 * $this->p * $this->p + $this->K3 * pow($this->p, 4.0);
        if($latitude < 0.0)
            $northing = 10000000.0 + $northing;
        return $northing;
    }

    protected function getEasting() {
        return 500000.0 + ($this->K4 * $this->p + $this->K5 * pow($this->p, 3.0));
    }
} // LatLon2UTM

class Digraphs {
    private $digraph1;
    private $digraph2;
    private $digraph1Array = array("A", "B", "C", "D", "E", "F", "G", "H", "J", "K", "L", "M", "N", "P",
                                       "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z");
    private $digraph2Array = array("V", "A", "B", "C", "D", "E", "F", "G", "H", "J", "K", "L", "M", "N",
                                       "P", "Q", "R", "S", "T", "U", "V");

    function __construct() {
        for($i=1; $i<25; $i++)
            $this->digraph1[$i] = $this->digraph1Array[$i-1];
         //   $this->digraph1($i => $digraph1Array[$i-1]);

        for($i=0; $i<=20; $i++)
            $this->digraph2[$i] = $this->digraph2Array[$i];
        //    $this->digraph2($i => $digraph2Array[$i]);
    }

    public function getDigraph1Index($letter) {
       assert('!empty($letter) && isset($letter)', 'input variable letter is empty');
       for($i=0; $i<count($this->digraph1Array); $i++)
           if($this->digraph1Array[$i]==$letter)
               return $i+1;
       return -1;
    }

    public function getDigraph2Index($letter) {
        assert('!empty($letter) && isset($letter)', 'input variable letter is empty');
        for($i=0; $i<count($this->digraph2Array); $i++)
            if($this->digraph2Array[$i]==$letter)
                return $i;
        return -1;
    }

    public function getDigraph1($longZone, $easting) {
        assert('!empty($longZone) && isset($longZone)', 'input variable longZone is empty');
        assert('!empty($easting) && isset($easting)', 'input variable easting is empty');
        $a1 = $longZone;
        $a2 = 8.0 * (($a1 - 1.0) % 3.0) + 1.0;
        $a3 = $easting;
        $a4 = $a2 + ((int)($a3/100000.0)) - 1.0;
        return (string)($this->digraph1[(int)(floor($a4))]);
    }

    public function getDigraph2($longZone, $northing) {
        assert('!empty($longZone) && isset($longZone)', 'input variable longZone is empty');
        assert('!empty($northing) && isset($northing)', 'input variable northing is empty');
        $a1 = $longZone;
        $a2 = 1.0 + 5.0 * (($a1-1.0)%2.0);
        $a3 = $northing;
        $a4 = ($a2 + ((int)($a3/100000.0)));
        $a4 = ($a2 + ((int)($a3 / 100000.0)))%20;
        $a4 = floor($a4);
        if($a4 < 0)
            $a4 = $a4 + 19.0;
        return (string)($this->digraph2[(int)(floor($a4))]);
    }
}

class MGRUTM2LatLon extends UTM2LatLon {
    public function convertMGRUTMToLatLong($mgrutm) {
       assert('!empty($mgrutm) && isset($mgrutm)', 'input variable mgrutm is empty');
       $zone = intval(substr($mgrutm, 0, 2));
       $latZone = substr($mgrutm, 2, 1);
       $digraph1 = substr($mgrutm, 3, 1);
       $digraph2 = substr($mgrutm, 4, 1);
       $this->easting = doubleval(substr($mgrutm, 5, 5));
       $this->northing = doubleval(substr($mgrutm, 10, 5));

       $lz = new LatZones();
       $latZoneDegree = $lz->getLatZoneDegree($latZone);

       $a1 = $latZoneDegree * 40000000.0 / 360.0;
       $a2 = 2000000 * floor($a1 / 2000000.0);

       $digraphs = new Digraphs();

       $digraph2Index = $digraphs->getDigraph2Index($digraph2);

       $startindexEquator = 1;
       if((1+$zone%2)==1)
           $startindexEquator = 6.0;

       $a3 = $a2 + ($digraph2Index - $startindexEquator) * 100000;
       if($a3<=0)
           $a3 = 10000000 + $a3;
       $this->northing = $a3 + $this->northing;

       $this->zoneCM = -183 + 6 * $zone;
       $digraph1Index = $digraphs->getDigraph1Index($digraph1);
       $a5 = 1 + $zone % 3;
       $a6 = array(16, 0, 8);
       $a7 = 100000 * ($digraph1Index - $a6[$a5 - 1]);
       $this->easting = $this->easting + $a7;
       $this->setVariables();
       $latitude = 180.0 * ($this->phi1 - $this->fact1 * ($this->fact2 + $this->fact3 + $this->fact4)) / pi();
       if($latZoneDegree < 0) 
            $latitude = 90.0 - $latitude;                                                           
       $d = $this->_a2 * 180.0 / pi();
       $longitude = $this->zoneCM - $d;
       if($this->getHemisphere($latZone)=="S") 
           $latitude = -$latitude;
       return array($latitude, $longitude);
    }
}

class LatLon2MGRUTM extends LatLon2UTM {

    public function validate($latitude, $longitude) {
        assert('!empty($latitude) && isset($latitude)', 'input variable latitude is empty');
        assert('!empty($longitude) && isset($longitude)', 'input variable longitude is empty');
        assert('$latitude >= -90.0 || $latitude <= 90.0 || $longitude >= -180.0 || $longitude < 180.0)', 
               'Legal ranges: latitude [-90,90], longitude [-180, 180)');
    }

    public function convertLatLonToMMRUTM( $latitude, $longitude) {
        assert('!empty($latitude) && isset($latitutde)', 'input variable lat is empty');
        assert('!empty($longitude) && isset($longitude)', 'input variable lon is empty');
    
        validate($latitude, $longitude);
        $this->setVariables($latitude, $longitude);

        $longZone = $this->getLongZone($longitude);
        $latZones = new LatZones();
        $latZone = $latZones.getLatZone($latitude);

        $_easting = $this->getEasting();
        $_northing = $this->getNorthing();

        $digraphs = new Digraphs();
        $digraph1 = $digraphs.getDigraph1(intval($longZone), $_easting);
        $digraph2 = $digraphs.getDigraph2(intval($longZone), $_northing);

        $easting = (string)((int)_easting);
        if(strlen($easting) < 5)
            $easting = "00000" . $easting;
        $easting = substr($easting, strlen($easting)-5);

        $northing = (string)((int)_northing);
        if(strlen($northing) < 5)
            $northing = "0000" . $northing;
        $northing = substr($northing, strlen($northing)-5);

        $mgrUTM = $longZone . $latZone . $digraph1 . $digraph2 . $easting . $northing;
        return $mgrUTM;
    }
}

/**
 * Main Coordinate Conversion Class
 */
class CoordConversion {
    public function utmToLatLon($utm) {
        assert('!empty($utm) && isset($utm)', 'input variable is empty');
        $c = new UTM2LatLon();
        return $c->convertUTMToLatLong($utm);
    }

    public function latlon2utm($latitude, $longitude) {
        assert('!empty($latitude) && isset($latitude)', 'input variable latitude is empty');
        assert('!empty($longitude) && isset($longitude)', 'input variable longitude is empty');
        $c = new LatLon2UTM();
        return $c->convertLatLonToUTM($latitude, $longitude);
    }

    public function latLon2MGRUTM($latitude, $longitude) {
        assert('!empty($latitude) && isset($latitude)', 'input variable latitude is empty');
        assert('!empty($longitude) && isset($longitude)', 'input variable longitude is empty');
        $c = new LatLon2MGRUTM();
        return $c->convertLatLonToMGRUTM($latitude, $longitude);
    }

    public function mgrutm2LatLon($MGRUTM) {
        assert('!empty($MGRUTM) && isset($MGRUTM)', 'input variable MGRUTM is empty');
        $c = new MGRUTM2LatLon();
        return $c->convertMGRUTMToLatLong($MGRUTM);
    }

}
