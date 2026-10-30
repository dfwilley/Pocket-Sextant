/* Pocket Sextant — navigation math core.
 * Sun ephemeris (Meeus low-precision, ~0.01 deg), sight corrections,
 * sight reduction (Hc/Zn), and a least-squares celestial fix with
 * running-fix support. Angles in degrees unless noted; intercepts in
 * arcminutes (= nautical miles).
 */
(function (root) {
  'use strict';
  var D2R = Math.PI / 180, R2D = 180 / Math.PI;
  function sin(d) { return Math.sin(d * D2R); }
  function cos(d) { return Math.cos(d * D2R); }
  function tan(d) { return Math.tan(d * D2R); }
  function asin(x) { return Math.asin(Math.max(-1, Math.min(1, x))) * R2D; }
  function atan2(y, x) { return Math.atan2(y, x) * R2D; }
  function norm360(a) { a = a % 360; return a < 0 ? a + 360 : a; }
  function norm180(a) { a = norm360(a); return a > 180 ? a - 360 : a; }

  /* ---------- Sun ephemeris ---------- */
  function julianDate(ms) { return ms / 86400000 + 2440587.5; }

  function sun(ms) {
    var JD = julianDate(ms), T = (JD - 2451545.0) / 36525;
    var L0 = norm360(280.46646 + T * (36000.76983 + 0.0003032 * T));
    var M = norm360(357.52911 + T * (35999.05029 - 0.0001537 * T));
    var e = 0.016708634 - T * (0.000042037 + 0.0000001267 * T);
    var C = (1.914602 - T * (0.004817 + 0.000014 * T)) * sin(M) +
            (0.019993 - 0.000101 * T) * sin(2 * M) + 0.000289 * sin(3 * M);
    var trueLong = L0 + C, nu = M + C;
    var R = 1.000001018 * (1 - e * e) / (1 + e * cos(nu));
    var omega = 125.04 - 1934.136 * T;
    var lambda = trueLong - 0.00569 - 0.00478 * sin(omega);
    var eps0 = 23 + (26 + (21.448 - T * (46.815 + T * (0.00059 - T * 0.001813))) / 60) / 60;
    var eps = eps0 + 0.00256 * cos(omega);
    var ra = norm360(atan2(cos(eps) * sin(lambda), cos(lambda)));
    var dec = asin(sin(eps) * sin(lambda));
    var gmst = norm360(280.46061837 + 360.98564736629 * (JD - 2451545.0) +
                       T * T * (0.000387933 - T / 38710000));
    var gast = gmst + (-0.00478 * sin(omega)) * cos(eps); // equation of the equinoxes
    var gha = norm360(gast - ra);
    var sd = 0.2666 / R;                 // semi-diameter, deg
    var hp = (8.794 / 3600) / R;         // horizontal parallax, deg
    // Equation of time: apparent minus mean solar time, in minutes
    var utDeg = ((ms % 86400000) + 86400000) % 86400000 / 3600000 * 15;
    var eot = norm180(gha - (utDeg - 180)) * 4;
    return { gha: gha, dec: dec, sd: sd, hp: hp, R: R, eotMin: eot };
  }

  /* ---------- Sight corrections ---------- */
  // Refraction in degrees for apparent altitude h (deg). Bennett (1982).
  function refraction(h, tempF, pressInHg) {
    if (h < -1) h = -1;
    var Rmin = 1 / tan(h + 7.31 / (h + 4.4));
    var tC = (tempF == null ? 50 : (tempF - 32) * 5 / 9);
    var pMb = (pressInHg == null ? 29.92 : pressInHg) * 33.8639;
    Rmin *= (pMb / 1010) * (283 / (273 + tC));
    return Math.max(0, Rmin) / 60;
  }
  function dip(heFt) { return heFt > 0 ? 0.97 * Math.sqrt(heFt) / 60 : 0; }

  /* sight: { t, hs (deg), method:'phone'|'marine', limb:'center'|'lower'|'upper',
              heFt, ic (deg), tempF, pressInHg }  */
  function observedAltitude(s) {
    var S = sun(s.t);
    var ha = s.hs + (s.ic || 0);
    if (s.method === 'marine') ha -= dip(s.heFt || 0);
    var R = refraction(ha, s.tempF, s.pressInHg);
    var pa = S.hp * cos(ha);
    var sdc = 0;
    if (s.limb === 'lower') sdc = S.sd; else if (s.limb === 'upper') sdc = -S.sd;
    return { ho: ha - R + pa + sdc, ha: ha, refr: R, pa: pa, sd: sdc,
             dip: s.method === 'marine' ? dip(s.heFt || 0) : 0, sun: S };
  }

  /* ---------- Sight reduction ---------- */
  function hcZn(lat, lon, gha, dec) {
    var lha = norm360(gha + lon);
    var hc = asin(sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(lha));
    var zn = norm360(atan2(-cos(dec) * sin(lha), cos(lat) * sin(dec) - sin(lat) * cos(dec) * cos(lha)));
    return { hc: hc, zn: zn, lha: lha };
  }

  // Move a position back along course/speed by dtHours (dead reckoning, plane sailing).
  function backtrack(lat, lon, courseDeg, speedKt, dtHours) {
    if (!speedKt || !dtHours) return { lat: lat, lon: lon };
    var d = speedKt * dtHours; // nm
    var la = lat - d * cos(courseDeg) / 60;
    var lo = lon - d * sin(courseDeg) / (60 * Math.max(0.01, cos((lat + la) / 2)));
    return { lat: la, lon: norm180(lo) };
  }
  function distanceNm(a, b) {
    var c = sin(a.lat) * sin(b.lat) + cos(a.lat) * cos(b.lat) * cos(b.lon - a.lon);
    return Math.acos(Math.max(-1, Math.min(1, c))) * R2D * 60;
  }
  function bearing(a, b) {
    return norm360(atan2(sin(b.lon - a.lon) * cos(b.lat),
      cos(a.lat) * sin(b.lat) - sin(a.lat) * cos(b.lat) * cos(b.lon - a.lon)));
  }

  /* ---------- Fix ---------- */
  // sights: array with .t and .ho (deg). opts: {course, speed, dr:{lat,lon}|null, priorSigmaMin}
  function solveFix(sights, opts) {
    opts = opts || {};
    var n = sights.length;
    if (n < 2) return { error: 'Need at least two sights for a fix.' };
    var course = opts.course || 0, speed = opts.speed || 0;
    var tf = Math.max.apply(null, sights.map(function (s) { return s.t; }));
    var pre = sights.map(function (s) {
      var S = sun(s.t);
      return { gha: S.gha, dec: S.dec, ho: s.ho, dt: (tf - s.t) / 3600000 };
    });
    function resid(lat, lon) {
      var r = new Array(n);
      for (var i = 0; i < n; i++) {
        var p = pre[i], q = backtrack(lat, lon, course, speed, p.dt);
        r[i] = (p.ho - hcZn(q.lat, q.lon, p.gha, p.dec).hc) * 60;
      }
      return r;
    }
    function cost(lat, lon) {
      var r = resid(lat, lon), c = 0;
      for (var i = 0; i < n; i++) c += r[i] * r[i];
      return c;
    }
    // Coarse global grid search for local minima
    var step = 1.5, lats = [], lons = [];
    for (var a = -84; a <= 84.001; a += step) lats.push(a);
    for (var o = -180; o < 180; o += step) lons.push(o);
    var grid = lats.map(function (la) { return lons.map(function (lo) { return cost(la, lo); }); });
    var cands = [];
    for (var i = 0; i < lats.length; i++) {
      for (var j = 0; j < lons.length; j++) {
        var c0 = grid[i][j], isMin = true;
        for (var di = -1; di <= 1 && isMin; di++) for (var dj = -1; dj <= 1; dj++) {
          if (!di && !dj) continue;
          var ii = i + di, jj = (j + dj + lons.length) % lons.length;
          if (ii < 0 || ii >= lats.length) continue;
          if (grid[ii][jj] < c0) { isMin = false; break; }
        }
        if (isMin) cands.push({ lat: lats[i], lon: lons[j], c: c0 });
      }
    }
    cands.sort(function (x, y) { return x.c - y.c; });
    cands = cands.slice(0, 12);
    if (opts.dr) cands.push({ lat: opts.dr.lat, lon: opts.dr.lon, c: cost(opts.dr.lat, opts.dr.lon) });

    function refine(p) {
      var lat = p.lat, lon = p.lon, lam = 1e-3, c = cost(lat, lon);
      for (var it = 0; it < 60; it++) {
        var r = resid(lat, lon), h = 1e-4;
        var r1 = resid(lat + h, lon), r2 = resid(lat, lon + h);
        var a11 = 0, a12 = 0, a22 = 0, b1 = 0, b2 = 0;
        for (var k = 0; k < n; k++) {
          var j1 = (r1[k] - r[k]) / h, j2 = (r2[k] - r[k]) / h;
          a11 += j1 * j1; a12 += j1 * j2; a22 += j2 * j2; b1 += j1 * r[k]; b2 += j2 * r[k];
        }
        var improved = false;
        for (var tries = 0; tries < 8; tries++) {
          var m11 = a11 * (1 + lam), m22 = a22 * (1 + lam), det = m11 * m22 - a12 * a12;
          if (Math.abs(det) < 1e-30) { lam *= 10; continue; }
          var d1 = -(m22 * b1 - a12 * b2) / det, d2 = -(m11 * b2 - a12 * b1) / det;
          // cap step size to 5 degrees
          var mag = Math.sqrt(d1 * d1 + d2 * d2); if (mag > 5) { d1 *= 5 / mag; d2 *= 5 / mag; }
          var nl = Math.max(-89.5, Math.min(89.5, lat + d1)), no = norm180(lon + d2);
          var nc = cost(nl, no);
          if (nc < c) { lat = nl; lon = no; c = nc; lam = Math.max(1e-7, lam / 10); improved = true; break; }
          lam *= 10;
        }
        if (!improved || Math.abs(d1) + Math.abs(d2) < 1e-9) break;
      }
      return { lat: lat, lon: lon, c: c };
    }
    var sols = [];
    cands.forEach(function (p) {
      var s = refine(p);
      for (var k = 0; k < sols.length; k++) if (distanceNm(sols[k], s) < 3) { if (s.c < sols[k].c) sols[k] = s; return; }
      sols.push(s);
    });
    sols.sort(function (x, y) { return x.c - y.c; });

    var prior = opts.priorSigmaMin || 15;
    var best = sols[0];
    var tol = best.c + Math.max(3 * best.c, n * prior * prior * 0.5);
    var plaus = sols.filter(function (s) { return s.c <= tol; });
    var chosen = best;
    if (opts.dr && plaus.length > 1) {
      chosen = plaus.slice().sort(function (x, y) { return distanceNm(x, opts.dr) - distanceNm(y, opts.dr); })[0];
    }
    var alternates = plaus.filter(function (s) { return s !== chosen; });

    // Residuals, azimuths, intercepts at chosen fix
    var lines = pre.map(function (p, i) {
      var q = backtrack(chosen.lat, chosen.lon, course, speed, p.dt);
      var r = hcZn(q.lat, q.lon, p.gha, p.dec);
      return { zn: r.zn, hc: r.hc, intercept: (p.ho - r.hc) * 60 };
    });
    // Uncertainty (nm): Jacobian in nm units
    var rr = resid(chosen.lat, chosen.lon), hh = 1e-4;
    var rLat = resid(chosen.lat + hh, chosen.lon), rLon = resid(chosen.lat, chosen.lon + hh);
    var s11 = 0, s12 = 0, s22 = 0, ssr = 0;
    var kx = 60 * Math.max(0.01, cos(chosen.lat)), ky = 60;
    for (var k = 0; k < n; k++) {
      var jy = (rLat[k] - rr[k]) / hh / ky, jx = (rLon[k] - rr[k]) / hh / kx;
      s11 += jx * jx; s12 += jx * jy; s22 += jy * jy; ssr += rr[k] * rr[k];
    }
    var sigma = n > 2 ? Math.max(Math.sqrt(ssr / (n - 2)), prior * 0.25) : prior;
    var detS = s11 * s22 - s12 * s12, radius95 = null;
    if (detS > 1e-12) {
      var c11 = s22 / detS, c22 = s11 / detS, c12 = -s12 / detS;
      var tr = c11 + c22, ev = tr / 2 + Math.sqrt(Math.max(0, tr * tr / 4 - (c11 * c22 - c12 * c12)));
      radius95 = 2.45 * sigma * Math.sqrt(ev);
    }
    // Azimuth spread (lines are undirected: compare mod 180)
    var spread = 0;
    for (var a1 = 0; a1 < n; a1++) for (var a2 = a1 + 1; a2 < n; a2++) {
      var dz = Math.abs(norm180(lines[a1].zn - lines[a2].zn)) % 180;
      dz = Math.min(dz, 180 - dz); if (dz > spread) spread = dz;
    }
    return {
      lat: chosen.lat, lon: chosen.lon, fixTime: tf, lines: lines,
      rmsMin: Math.sqrt(ssr / n), sigmaMin: sigma, radius95: radius95,
      azSpread: spread, alternates: alternates.map(function (s) { return { lat: s.lat, lon: s.lon }; })
    };
  }

  /* ---------- Local apparent noon ---------- */
  function localNoon(lonDeg, dayMs) {
    // dayMs: any time on the UTC date of interest; returns ms of meridian passage
    var day0 = Math.floor(dayMs / 86400000) * 86400000;
    var t = day0 + (12 - lonDeg / 15) * 3600000;
    for (var i = 0; i < 4; i++) {
      var g = sun(t).gha, lha = norm180(g + lonDeg);
      t -= lha / 360.9856 * 86400000;
    }
    return t;
  }

  /* ---------- Formatting ---------- */
  function fmtDM(deg, pos, neg, degWidth) {
    var h = deg < 0 ? neg : pos, a = Math.abs(deg);
    var d = Math.floor(a), m = (a - d) * 60;
    if (m >= 59.95) { d += 1; m = 0; }
    var ds = String(d); while (ds.length < (degWidth || 0)) ds = '0' + ds;
    return h + ' ' + ds + '° ' + m.toFixed(1).padStart(4, '0') + '′';
  }
  function fmtAlt(deg) {
    var s = deg < 0 ? '−' : '', a = Math.abs(deg), d = Math.floor(a), m = (a - d) * 60;
    if (m >= 59.95) { d += 1; m = 0; }
    return s + d + '° ' + m.toFixed(1).padStart(4, '0') + '′';
  }
  // Parse "32 13.4 N", "N32°13.4'", "-110.9", "110 55.2 W"
  function parseCoord(str, isLat) {
    if (str == null) return null;
    var s = String(str).trim().toUpperCase();
    if (!s) return null;
    var neg = /[SW]/.test(s) || /^-/.test(s);
    var nums = s.replace(/[NSEW°'′"]/g, ' ').replace(/^-/, '').trim().split(/[\s,]+/).filter(Boolean).map(Number);
    if (!nums.length || nums.some(isNaN)) return null;
    var v = nums[0] + (nums[1] || 0) / 60 + (nums[2] || 0) / 3600;
    if (neg) v = -v;
    if (isLat ? Math.abs(v) > 90 : Math.abs(v) > 180) return null;
    return v;
  }

  var api = {
    sun: sun, refraction: refraction, dip: dip, observedAltitude: observedAltitude,
    hcZn: hcZn, backtrack: backtrack, distanceNm: distanceNm, bearing: bearing,
    solveFix: solveFix, localNoon: localNoon, fmtDM: fmtDM, fmtAlt: fmtAlt,
    parseCoord: parseCoord, norm180: norm180, norm360: norm360
  };
  root.SextantCore = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof self !== 'undefined' ? self : this);
