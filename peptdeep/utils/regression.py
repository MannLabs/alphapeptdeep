import pandas as pd
import numpy as np


def regional_sampling(
    psm_df: pd.DataFrame,
    target: str = "rt_norm",
    n_train: int = 1000,
    return_test_df: bool = False,
    random_state: int = 1337,
) -> pd.DataFrame:
    """Divide `psm_df` into 10 bins and sample training values
    from each bins for model fine-tuning. The values are defined in the `target`
    column (`rt_norm` or `ccs`).

    Parameters
    ----------
    psm_df : pd.DataFrame
        Dataframe of PSMs.

    target : str, optional
        Target columns to sample.
        Defaults to 'rt_norm'.

    n_train : int, optional
        The number of training PSMs to sample.
        Defaults to 1000.

    return_test_df : bool, optional
        If also return `test_df`.
        `test_df` contains the PSMs that are not sampled.
        Defaults to False.

    random_state : int
        `random_state` in `df.sample()`.

    Returns
    -------
    pd.DataFrame or tuple of pd.DataFrame
        The sampled training PSMs (DataFrame)
        Additional [pd.DataFrame] is returned if `return_test_df==True` for PSMs not sampled.
    """
    x = np.arange(0, 11) / 10 * psm_df[target].max()
    sub_n = n_train // (len(x) - 1)
    df_list = []
    for i in range(len(x) - 1):
        _df = psm_df[(psm_df[target] >= x[i]) & (psm_df[target] < x[i + 1])]
        if len(_df) == 0:
            pass
        elif len(_df) // 2 < sub_n:
            df_list.append(
                _df.sample(len(_df) // 2, replace=False, random_state=random_state)
            )
        else:
            df_list.append(_df.sample(sub_n, replace=False, random_state=random_state))
    if return_test_df:
        if len(df_list) == 0:
            return pd.DataFrame(), pd.DataFrame()
        train_df = pd.concat(df_list)
        test_df = psm_df.drop(train_df.index)
        return train_df, test_df
    else:
        if len(df_list) == 0:
            return pd.DataFrame()
        return pd.concat(df_list)


# legacy
uniform_sampling = regional_sampling


def linear_regression(x, y):
    coeffs = np.polyfit(x, y, 1)
    w, b = coeffs.tolist()
    yhat = np.poly1d(coeffs)(x)
    ybar = np.sum(y) / len(y)
    ssreg = np.sum((yhat - ybar) ** 2)
    sstot = np.sum((y - ybar) ** 2)
    R_square = ssreg / sstot
    return dict(
        R_square=[R_square],
        R=[np.sqrt(R_square)],
        slope=[w],
        intercept=[b],
    )


def evaluate_linear_regression(
    df: pd.DataFrame, x="rt_pred", y="rt_norm", ci=95, n_sample=10000000
):
    if len(df) > n_sample:
        df = df.sample(n_sample, replace=False)

    regs = linear_regression(df[x].values, df[y].values)
    regs["test_num"] = len(df)

    return pd.DataFrame(regs)


def evaluate_linear_regression_plot(
    df: pd.DataFrame, x="rt_pred", y="rt_norm", ci=95, n_sample=100000
):
    import seaborn as sns

    if len(df) > n_sample:
        df = df.sample(n_sample)
    alpha = 0.05
    if len(df) < 5000:
        alpha = 1
    elif len(df) < 50000:
        alpha = 5000.0 / len(df)
    return sns.regplot(
        data=df,
        x=x,
        y=y,
        color="r",
        ci=ci,
        scatter_kws={"s": 0.05, "alpha": alpha, "color": "b"},
    )
