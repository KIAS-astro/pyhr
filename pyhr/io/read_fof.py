import os
import pathlib
import struct
import numpy as np
import pandas as pd

class ReadFoF(object):
    def __init__(self, basedir='/ramses2/jaehyun/HR5/FoF_Data/', num=None, hid=None):
        self.basedir = pathlib.Path(basedir)
        self.fname_base = dict()
        self.fname_base['list'] = r'FoF.{0:s}/GALCATALOG.LIST.{0:s}'
        self.fname_base['data'] = r'FoF.{0:s}/GALFIND.DATA.{0:s}'

        self.num = num
        self.hid = hid

        # fname = self.basedir / fname_base['list']
        self.dtype, self.size = self._get_dtype_and_size()
        self.get_info(self.num)

    def get_info(self, num=None):
        if num is None:
            num = self.num
        elif num != self.num:
            self.num = num

        fname = self.basedir / self.fname_base['list'].format(str(num).zfill(5))
        print('Read nhalo')
        self.nhalo, self.nsub = self.read_catalogue(fname)
        print('Read halos')
        halos, subhalos, offsets = self.read_catalogue(fname, self.nhalo, self.nsub, nhalo_only=False)
        print('Converting to DataFrame')
        # pandas DataFrame containing halo and subhalo info
        df = pd.DataFrame(halos)
        dfs = pd.DataFrame(subhalos)
        df['offsets'] = offsets
        df['size_all'] = df['nDM']*self.size['DM'] + df['nstar']*self.size['star'] +\
                         df['nBH']*self.size['BH'] + df['ngas']*self.size['gas']

        # Byte offset
        df['offset2'] = (df['size_all'].cumsum()).shift(1, fill_value=0)
        # Starting subhalo id
        df['sid_start'] = df['nsub'].cumsum().shift(1, fill_value=0)
        columns_to_move = ['nsub', 'sid_start', 'offset2', 'size_all'] # + [f'size_{k}' for k in ['star', 'DM', 'gas', 'BH']]
        columns =  columns_to_move + [col for col in df.columns if col not in columns_to_move]
        df = df[columns]

        self.df = df
        self.dfs = dfs

        return df, dfs

    # Read all particles/grids of subhalos and unbound components of target FoF halos.
    def get_halo(self, hid):
        fname = self.basedir / self.fname_base['data'].format(str(self.num).zfill(5))

        df = self.df.loc[hid].copy()
        dfs = self.dfs

        dat = dict()
        nelem = dict()
        indices = dict()
        for i in df.index:
            dat[i] = dict()
            nelem[i] = dict()
            indices[i] = dict()
            for k in ['DM', 'star', 'BH', 'gas']:
                dat[i][k] = np.zeros(df.loc[i, f'n{k}'], dtype=self.dtype[k])

            sid0 = df.loc[i, f'sid_start']
            nsub = df.loc[i, f'nsub']
            for k in ['DM', 'star', 'BH', 'gas']:
                nelem[i][k] = []
            for isub in range(nsub):
                for k in ['DM', 'star', 'BH', 'gas']:
                    nelem[i][k].append(dfs.loc[sid0 + isub, f'n{k}'])

            for k in ['DM', 'star', 'BH', 'gas']:
                nelem[i][k] = np.array(nelem[i][k])
                indices[i][k] = np.insert(np.cumsum(nelem[i][k]), 0, 0)

            with open(fname, 'rb') as fp:
                fp.seek(df.loc[i, 'offset2'], os.SEEK_SET)
                for k in ['DM', 'star', 'BH', 'gas']:
                    idx = indices[i][k]
                    for isub in range(nsub):
                        count = nelem[i][k][isub]
                        dat[i][k][idx[isub]:idx[isub+1]] = np.frombuffer(
                            fp.read(count*self.size[k]), dtype=self.dtype[k],
                            count=count)

        return dat, nelem, indices

    def read_catalogue(self, fname, nhalo=None, nsub=None, nhalo_only=False):
        if nhalo is None:
            nhalo_only = True

        if not nhalo_only:
            # Halo info
            halos = np.empty(nhalo, dtype=self.dtype['halo'])
            subhalos = np.empty(nsub, dtype=self.dtype['subhalo'])
            # Starting file position of halos
            offsets = np.empty(nhalo, dtype=np.int64)

        with open(fname, 'rb') as fp:
            # Find EOF
            fp.seek(0, os.SEEK_END)
            eof = fp.tell()
            # Go to the start of the stream
            fp.seek(0, os.SEEK_SET)

            nhalo = 0
            nsub = 0
            while fp.tell() < eof:
                offset = fp.tell()
                halo = np.frombuffer(fp.read(self.size['halo']), dtype=self.dtype['halo'])
                nsub_ = halo['nsub'].item()
                # Record halo info and start position of subhalos
                if not nhalo_only:
                    subhalo_ = np.frombuffer(fp.read(nsub_*self.size['subhalo']), dtype=self.dtype['subhalo'],
                                             count=nsub_)
                    halos[nhalo] = halo
                    subhalos[nsub:nsub+nsub_] = subhalo_
                    offsets[nhalo] = offset
                else:
                    # Jump to the next halo
                    fp.seek(nsub_*self.size['subhalo'], os.SEEK_CUR)

                nhalo += 1
                nsub += nsub_

        if not nhalo_only:
            return halos, subhalos, offsets
        else:
            return nhalo, nsub

    def _get_dtype_and_size(self):
        """
        For more details, see
        https://astro.kias.re.kr/Horizon-Run5/index.php/FoF_halo/substructure_finding_data
        """

        dtype = dict()
        # Define data type
        # Halo and subhalo (GALCATALOG.LIST)
        dtype['halo'] = np.dtype([('nsub', 'i4'), ('nDM', 'i4'), ('nstar', 'i4'),
                                  ('nBH', 'i4'), ('ngas', 'i4'), ('npall', 'i4'),
                                  ('mtot', 'f8'), ('mdm', 'f8'), ('mgas', 'f8'),
                                  ('mbh', 'f8'), ('mstar', 'f8'),
                                  ('x1','f8'), ('x2','f8'), ('x3','f8'),
                                  ('v1','f8'), ('v2','f8'), ('v3','f8')])

        dtype['subhalo'] = np.dtype([('nDM', 'i4'), ('ngas', 'i4'), ('nBH', 'i4'),
                                     ('nstar', 'i4'), ('npall', 'i4'), ('dum', 'i4'),
                                     ('mtot', 'f8'), ('mdm', 'f8'), ('mgas', 'f8'),
                                     ('mbh', 'f8'), ('mstar', 'f8'),
                                     ('x1','f8'), ('x2','f8'), ('x3','f8'),
                                     ('v1','f8'), ('v2','f8'), ('v3','f8')])

        # Particle (DM, star, BH) and gas
        dtype['DM'] = np.dtype([('x1', 'f8'), ('x2', 'f8'), ('x3', 'f8'),
                                ('v1', 'f8'), ('v2', 'f8'), ('v3', 'f8'),
                                ('mass', 'f8'), ('dum0', 'f8'), ('tp', 'f8'),
                                ('zp', 'f8'), ('mass0', 'f8'), ('tpp', 'f8'),
                                ('indtab', 'f8'), ('id', 'i8'), ('potential', 'f8'),
                                ('level', 'i4'), ('dum1', 'i4')])

        dtype['star'] = dtype['DM']

        dtype['BH'] = np.dtype([('x1', 'f8'), ('x2', 'f8'), ('x3', 'f8'),
                                ('v1', 'f8'), ('v2', 'f8'), ('v3', 'f8'),
                                ('mass', 'f8'), ('tbirth', 'f8'),
                                ('angm1', 'f8'), ('angm2', 'f8'), ('angm3', 'f8'),
                                ('ang1', 'f8'), ('ang2', 'f8'), ('ang3', 'f8'),
                                ('dm_real', 'f8'), ('dm_Bondi', 'f8'), ('dm_Edd', 'f8'),
                                ('esave', 'f8'), ('emag', 'f8'), ('eps', 'f8'),
                                ('id', 'i4'), ('dum', 'i4')])

        dtype['gas'] = np.dtype([('x1', 'f8'), ('x2', 'f8'), ('x3', 'f8'), ('dx', 'f8'),
                                 ('v1', 'f8'), ('v2', 'f8'), ('v3', 'f8'), ('dum0', 'f4'),
                                 ('den', 'f8'), ('press', 'f4'), ('Z', 'f4'),
                                 ('fe', 'f4'), ('h', 'f4'), ('o', 'f4'),
                                 ('level', 'i4'), ('mass', 'f4'), ('dum1', 'f4'),
                                 ('id', 'i8'), ('potential', 'f8'),
                                 ('f1', 'f8'), ('f2', 'f8'), ('f3', 'f8')])

        size = {k: v.itemsize for k, v in dtype.items()}

        return dtype, size
