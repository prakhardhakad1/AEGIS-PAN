/**
 * AEGIS-PAN v5.0 Isomorphic JavaScript Competing-Risks & Royston-Parmar Engine (`web/js/engine.mjs`).
 * Evaluates the exact Royston-Parmar natural cubic spline basis v_1(log t), v_2(log t) and
 * Fine-Gray competing-risks cumulative incidence functions from `results/model_artifact.json`.
 */

export function roystonParmarBasisJS(logT, knots) {
  const [kMin, kInt, kMax] = knots;
  const lam = (kMax - kInt) / (kMax - kMin);
  const posInt = Math.max(0.0, logT - kInt);
  const posMin = Math.max(0.0, logT - kMin);
  const posMax = Math.max(0.0, logT - kMax);
  const v1 = logT;
  const v2 = (Math.pow(posInt, 3) - lam * Math.pow(posMin, 3) - (1.0 - lam) * Math.pow(posMax, 3)) / Math.pow(kMax - kMin, 2);
  return [v1, v2];
}

export function evalRoystonParmarH0JS(m, rpGamma, rpKnots) {
  if (m <= 0) return 0.0;
  const [v1, v2] = roystonParmarBasisJS(Math.log(m), rpKnots);
  const logH0 = rpGamma[0] + rpGamma[1] * v1 + rpGamma[2] * v2;
  return Math.exp(Math.max(-20.0, Math.min(10.0, logH0)));
}

export function scoreProfileFromArtifactJS(art, prof, kappa = 1.0) {
  const org = prof.organ;
  const ageMean = art.organ_age_means[org];
  const ob = art.organ_baselines[org];
  const b1 = art.cause1_beta;
  const b2 = art.cause2_beta;

  const x1 = [
    (prof.age - ageMean) / 10.0,
    prof.T >= 2 ? 1.0 : 0.0,
    prof.T >= 3 ? 1.0 : 0.0,
    prof.T >= 4 ? 1.0 : 0.0,
    prof.N >= 1 ? 1.0 : 0.0,
    prof.N >= 2 ? 1.0 : 0.0,
    prof.grade - 2.0,
    prof.adjuvant ? 1.0 : 0.0
  ];
  const x2 = [
    (prof.age - 60.0) / 10.0,
    prof.comorb - 1.0
  ];

  let eta1 = Math.log(kappa);
  let eta1NoTx = Math.log(kappa);
  for (let i = 0; i < b1.length; i++) {
    eta1 += b1[i] * x1[i];
    if (i < b1.length - 1) eta1NoTx += b1[i] * x1[i];
  }
  const eta2 = b2[0] * x2[0] + b2[1] * x2[1];

  const shr1 = Math.exp(eta1);
  const hr2 = Math.exp(eta2);

  const curve = art.horizons.map(m => {
    if (m === 0) {
      return { m: 0, os: 1.0, cif1: 0.0, cif2: 0.0, cif1NoTx: 0.0, allCauseMort: 0.0 };
    }
    const h1Rp = evalRoystonParmarH0JS(m, art.cause1_rp_gamma, art.cause1_rp_knots) * ob.h1_scale;
    const h2Rp = evalRoystonParmarH0JS(m, art.cause2_rp_gamma, art.cause2_rp_knots) * ob.h2_scale;

    const h1 = h1Rp * shr1;
    const h1NoTx = h1Rp * Math.exp(eta1NoTx);
    const h2 = h2Rp * hr2;
    const hTot = h1 + h2;

    const osVal = Math.exp(-hTot);
    const totD = 1.0 - osVal;
    const cif1 = totD * (h1 / Math.max(1e-12, hTot));
    const cif2 = totD * (h2 / Math.max(1e-12, hTot));

    const osNoTx = Math.exp(-(h1NoTx + h2));
    const cif1NoTx = (1.0 - osNoTx) * (h1NoTx / Math.max(1e-12, h1NoTx + h2));

    return {
      m,
      os: osVal,
      cif1,
      cif2,
      cif1NoTx,
      allCauseMort: totD
    };
  });

  return { shr1, hr2, curve };
}
