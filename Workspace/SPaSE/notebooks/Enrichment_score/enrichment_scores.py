import pandas as pd
import numpy as np

def score_abid_pval(df, program, nbins, nbackground, obs_cut = None):
  if obs_cut is None:
    obs_avg = df.mean(axis=0)
    n_items = int(np.round(len(obs_avg) / (nbins - 1)))
    obs_cut = (obs_avg.rank(method='min') // n_items).astype("int")


  program_mean = df.T.loc[program].mean(axis=0)

  background = []
  for i in range(nbackground):
      background_genes = []

      for cut in obs_cut.loc[program]:
          r_genes = set(obs_cut[obs_cut == cut].index) - set(program)
          
          r_genes = np.array(list(r_genes))
          np.random.shuffle(r_genes)
          # uses full r_genes if ctrl_size > len(r_genes)
          background_genes.append(r_genes[0])

      # get the mean expression of the background genes
      background.append(df.T.loc[background_genes].mean(axis=0))

  background = pd.DataFrame(background).T
  result = (program_mean < background[0]).astype(int)
  for i in range(1, len(background.columns)):
      result = result + (program_mean < background[i]).astype(int)

  pval = result / len(background.columns)

  return pval


def score_abid_enrichment(df, program, nbins, nbackground, obs_cut = None):
  if obs_cut is None:
    obs_avg = df.mean(axis=0)
    n_items = int(np.round(len(obs_avg) / (nbins - 1)))
    obs_cut = (obs_avg.rank(method='min') // n_items).astype("int")


  program_mean = df.T.loc[program].mean(axis=0)

  background = []
  for i in range(nbackground):
      background_genes = []

      for cut in obs_cut.loc[program]:
          r_genes = set(obs_cut[obs_cut == cut].index) - set(program)
          
          r_genes = np.array(list(r_genes))
          np.random.shuffle(r_genes)
          # uses full r_genes if ctrl_size > len(r_genes)
          background_genes.append(r_genes[0])

      # get the mean expression of the background genes
      background.append(df.T.loc[background_genes].mean(axis=0))

  background = pd.DataFrame(background).T
  result = program_mean - background.mean(axis=1)

  return result


def score_abid_zscore(df, program, nbins, nbackground, obs_cut = None):
  if obs_cut is None:
    obs_avg = df.mean(axis=0)
    n_items = int(np.round(len(obs_avg) / (nbins - 1)))
    obs_cut = (obs_avg.rank(method='min') // n_items).astype("int")


  program_df = df.T.loc[program].mean(axis=0)

  background = []
  for i in range(nbackground):
      background_genes = []

      for cut in obs_cut.loc[program]:
          r_genes = set(obs_cut[obs_cut == cut].index) - set(program)
          
          r_genes = np.array(list(r_genes))
          np.random.shuffle(r_genes)
          # uses full r_genes if ctrl_size > len(r_genes)
          background_genes.append(r_genes[0])

      # get the mean expression of the background genes
      background.append(df.T.loc[background_genes].T)

  background = pd.concat(background)
  background_mean = background.groupby(background.index).mean()
  background_std = background.groupby(background.index).std()

  z_scores = (program_df - background_mean)/(background_std + 1e-12)

  result = z_scores.sum(axis=1)

  return result


def score_abid_pval_and_enrichment(df, program, nbins, nbackground, obs_cut = None):
  if obs_cut is None:
    obs_avg = df.mean(axis=0)
    n_items = int(np.round(len(obs_avg) / (nbins - 1)))
    obs_cut = (obs_avg.rank(method='min') // n_items).astype("int")


  program_mean = df.T.loc[program].mean(axis=0)

  background = []
  for i in range(nbackground):
      background_genes = []

      for cut in obs_cut.loc[program]:
          r_genes = set(obs_cut[obs_cut == cut].index) - set(program)
          
          r_genes = np.array(list(r_genes))
          np.random.shuffle(r_genes)
          # uses full r_genes if ctrl_size > len(r_genes)
          background_genes.append(r_genes[0])

      # get the mean expression of the background genes
      background.append(df.T.loc[background_genes].mean(axis=0))

  background = pd.DataFrame(background).T
  result_enrichment = program_mean - background.mean(axis=1)
  
  result_pval = (program_mean < background[0]).astype(int)
  for i in range(1, len(background.columns)):
      result_pval = result_pval + (program_mean < background[i]).astype(int)

  pval = result_pval / len(background.columns)

  return pval, result_enrichment
