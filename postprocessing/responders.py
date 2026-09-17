"""Profile-level and exclusive responder classifications."""
import numpy as np
import pandas as pd

X1_CONCENTRATIONS_UM = {"flecainide": 1.5, "vernakalant": 10.0, "tertiapin-Q": 0.1}


def exclusive_responders(results: pd.DataFrame) -> pd.DataFrame:
    """Derive Figure 6 exclusive responders from the public long-format table."""
    required = {"population", "profile_id", "drug", "concentration_uM", "profile_responder"}
    missing = required.difference(results.columns)
    if missing:
        raise ValueError(f"missing required columns: {sorted(missing)}")
    selected = []
    for drug, concentration in X1_CONCENTRATIONS_UM.items():
        part = results.loc[
            (results["drug"] == drug)
            & np.isclose(pd.to_numeric(results["concentration_uM"], errors="coerce"), concentration)
        ].copy()
        part["response"] = pd.to_numeric(part["profile_responder"], errors="coerce")
        # Repeated S2 rows carry the same profile flag; max is also robust to sparse rows.
        part = part.groupby(["population", "profile_id"], as_index=False)["response"].max()
        part["drug"] = drug
        selected.append(part)
    wide = pd.concat(selected, ignore_index=True).pivot(
        index=["population", "profile_id"], columns="drug", values="response"
    ).reset_index()
    for drug in X1_CONCENTRATIONS_UM:
        if drug not in wide:
            wide[drug] = np.nan
    drug_columns = list(X1_CONCENTRATIONS_UM)
    complete = wide[drug_columns].notna().all(axis=1)
    positive_count = wide[drug_columns].fillna(0).eq(1).sum(axis=1)
    wide["exclusive_drug"] = pd.NA
    for drug in drug_columns:
        mask = complete & positive_count.eq(1) & wide[drug].eq(1)
        wide.loc[mask, "exclusive_drug"] = drug
    wide["is_exclusive_responder"] = wide["exclusive_drug"].notna()
    return wide.sort_values(["population", "profile_id"], kind="stable").reset_index(drop=True)