import os
import functools

import os.path as osp
import pandas as pd
from inspect import getcallargs

def pickle_galfind_list(func):

    @functools.wraps(func)
    def wrapper(cls, *args, **kwargs):
        # Convert positional args to keyword args
        call_args = getcallargs(func, cls, *args, **kwargs)
        call_args.pop('self')
        kwargs = call_args
        name = func.__name__

        if 'num' not in kwargs:
            kwargs['num'] = cls.num
        else:
            cls.num = kwargs['num']

        num = kwargs['num']

        if 'prefix' in kwargs:
            prefix = kwargs['prefix']
        else:
            prefix = None

        if 'savdir' in kwargs:
            savdir = kwargs['savdir']
            if savdir is None:
                savdir = cls.savdir
        else:
            savdir = cls.savdir

        if prefix is not None:
            savdir = osp.join(savdir, prefix)

        if 'force_override' in kwargs:
            force_override = kwargs['force_override']
        else:
            force_override = False

        # Create savdir if it doesn't exist
        if not osp.exists(savdir):
            try:
                os.makedirs(savdir)
                force_override = True
            except (IOError, PermissionError) as e:
                cls.logger.warning('Could not make directory')

        if prefix is None:
            fpkl1 = osp.join(savdir, 'FoF.GALCATALOG.{0:05d}.halo.p'.\
                             format(num))
            fpkl2 = osp.join(savdir, 'FoF.GALCATALOG.{0:05d}.subhalo.p'.\
                             format(num))
        else:
            fpkl1 = osp.join(savdir, 'FoF.GALCATALOG.{0:05d}.{1:s}.halo.p'.\
                             format(num, prefix))
            fpkl2 = osp.join(savdir, 'FoF.GALCATALOG.{0:05d}.{1:s}.subhalo.p'.\
                             format(num, prefix))

        # Check if the original catalog file is updated
        if (osp.exists(fpkl1) and osp.exists(fpkl2)) and not force_override:
            cls.logger.info(f'[{name}]: Reading pickle..')
            cls.df = pd.read_pickle(fpkl1)
            cls.dfs = pd.read_pickle(fpkl2)
            return cls.df, cls.dfs
        else:
            cls.logger.info(f'[{name}]: Reading original file..')
            # Call function
            df, dfs = func(cls, *args, **kwargs)
            try:
                df.to_pickle(fpkl1)
                dfs.to_pickle(fpkl2)
            except (IOError, PermissionError) as e:
                cls.logger.warning(
                    f'[{name}]: Cannot not pickle to {0:s}.'.format(fpkl1))
            return df, dfs

    return wrapper
