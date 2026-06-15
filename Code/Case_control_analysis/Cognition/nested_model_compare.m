function out = nested_model_compare(mdl_small, mdl_large)
    RSS_small = mdl_small.SSE;
    RSS_large = mdl_large.SSE;
    df_small = mdl_small.DFE;
    df_large = mdl_large.DFE;

    df1 = df_small - df_large;
    df2 = df_large;

    F = ((RSS_small - RSS_large) / df1) / (RSS_large / df2);
    out.F = F;
    out.df1 = df1;
    out.df2 = df2;
    out.pValue = 1 - fcdf(F, df1, df2);
end



