# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2017  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Maria Fedorova
# LPPtiger is Dual-licensed
#     For academic and non-commercial use: `GPLv2 License` Please read more information by the following link:
#         [The GNU General Public License version 2] (https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
#     For commercial use:
#         please contact the SysMedOs_team by email.
# Please cite our publication in an appropriate form.
#
# For more info please contact:
#     SysMedOs_team: oxlpp@bbz.uni-leipzig.de
#     LPPtiger repository: https://bitbucket.org/SysMedOs/lpptiger
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#

import os
import shutil
import pandas as pd


class LogPageCreator(object):
    def __init__(self, output_folder, start_time, params):
        print(os.getcwd())
        self.output_folder = output_folder
        self.output_img_folder = output_folder + r'/LPPtiger_Results_Figures_%s' % start_time
        self.main_page = output_folder + r'/LPPtiger_Results_%s.html' % start_time
        self.logo = r'LPPtiger_Results_Figures_%s/LPPtiger.ico' % start_time
        _image_lst_page = r'LPPtiger_Results_Figures_%s/LPPtiger_Results_Figures_list.html' % start_time
        _params_lst_page = r'LPPtiger_Results_Figures_%s/LPPtiger_Params_list.html' % start_time
        _idx_lst_page = r'LPPtiger_Results_Figures_%s/LPPtiger_Identification_list.html' % start_time
        self.image_lst_page = self.output_img_folder + r'/LPPtiger_Results_Figures_list.html'
        self.params_lst_page = self.output_img_folder + r'/LPPtiger_Params_list.html'
        self.idx_lst_page = self.output_img_folder + r'/LPPtiger_Identification_list.html'

        self.lipid_type = params['lipid_type']
        hunter_folder = params['hunter_folder']

        if params['fast_isotope'] is True:
            isotope_score_mode = '(Fast mode)'
        else:
            isotope_score_mode = ''

        with open(self.main_page, 'w') as _m_page:
            m_info_lst = ['<html>\n', '<link rel="icon" href="', self.logo, '" type="image/x-icon"/>\n'
                                                                            '<title>LPPtiger_Results ', start_time,
                          '</title>\n<frameset cols="390,*">\n<frameset rows="430,*">\n',
                          '<frame src="', _params_lst_page, '" frameborder="0" >\n',
                          '<frame src="', _idx_lst_page, '" frameborder="0" >\n</frameset>\n',
                          '<frame src="', _image_lst_page, '"name ="results_frame">\n</frameset>\n</html>\n'
                          ]
            _m_page.write(''.join(m_info_lst))

        with open(self.image_lst_page, 'w') as _img_page:
            _img_page.write('''
                            <html>\n<body>\n<style type="text/css">\n
                            p {margin-left: 16px; text-decoration: none; font-family: sans-serif;}\n
                            h3 {font-size:20px; margin-left: 16px; text-decoration: none; font-family: sans-serif;}\n
                            body {font-family: sans-serif;}\n
                            table, th, td {font-size:14px;text-align: center; font-family: sans-serif;}\n
                            th{background-color:#0066B2;color:white; margin:center;}\n
                            tr:nth-child(odd){background-color: #B1D3EC;}\n
                            tr:nth-child(even){background-color: #7C94A5;}\n
                            a:link {text-decoration:none; color:black} a:hover{text-decoration:underline; color:black;}
                            a:visited {text-decoration:none; color:black;}\n</style>\n''')

        with open(self.params_lst_page, 'w') as _params_page:
            _params_page.write('''
                                <html>\n<body>\n<style type="text/css">\n
                                p {margin-left: 20px; text-decoration: none; font-family: sans-serif;}\n
                                body {font-family: sans-serif;}\n table{width:100%s;}\n
                                table, th, td {font-size:14px;text-align: center; font-family: sans-serif;}\n
                                th{background-color:#0066B2;color:white; margin:center;}\n
                                a:link {text-decoration:none} a:hover{text-decoration:underline }\n
                                ul {font-size:14px; width: 260px;}\n </style>\n
                                <h3><img src="LPPtiger.ico" height=36/>  LPPtiger</h3><font size="1">\n
                                <hr> <h3>Parameters:</h3>\n<ul>\n
                                <li>Start time: %s</li>\n<li>Mode: %s %s</li>\n
                                <li><i>m/z</i> range: %.1f - %.1f <i>m/z</i></li>\n<li>RT range: %.1f - %.1f min</li>\n
                                <li>MS1 Threshold: %i</li>\n<li>MS2 Threshold: %i</li>\n
                                <li>MS1 ppm: %i</li>\n<li>MS2 ppm: %i</li>\n
                                <li>LPPtiger score > %.1f </li>\n<li>Isotope score > %.1f %s</li>\n
                                </ul>\n<hr>\n<h3>Lipid identification list:</h3><font size="1">\n<table>\n<thead>\n
                                <tr style="text-align: center;">\n
                                <th>ID#</th>\n<th> MS1_obs_mz </th>\n<th>RT(min)</th>\n<th>Discrete</th>\n
                                <th>Score</th>\n
                                </tr>\n</thead>\n</table>\n</body>\n</html>\n
                                ''' % ('%', params['hunter_start_time'], self.lipid_type, params['charge_mode'],
                                       params['mz_start'], params['mz_end'], params['rt_start'], params['rt_end'],
                                       params['ms_th'], params['ms2_th'], params['ms_ppm'], params['ms2_ppm'],
                                       params['score_filter'], params['isotope_score_filter'], isotope_score_mode))

        with open(self.idx_lst_page, 'w') as _idx_page:
            _idx_page.write('''
                            <html>\n<body>\n<style type="text/css">\n
                            p {margin-left: 20px; text-decoration: none; font-family: sans-serif;}\n
                            body {background-color: #B1D3EC;font-family: sans-serif;}\n table{width:100%s;}\n
                            table, th, td {font-size:14px;text-align: center; font-family: sans-serif;}\n
                            th{background-color:#0066B2;color:white;margin:center;}\n
                            tr:nth-child(even){background-color: #7C94A5;}\n
                            a:link {text-decoration:none; color:black}a:hover{text-decoration:underline; color:black;}
                            a:visited {text-decoration:none; color:black;}\n
                            </style>\n<table>\n<tbody>
                            ''' % '%')
        try:
            shutil.copy('%s\LPPtiger.ico' % hunter_folder, self.output_img_folder)
        except IOError:
            pass

    def add_info(self, img_path, ident_idx, ident_info_se):

        print('try to add identification to report html')

        ms1_pr_mz = ident_info_se['MS1_obs_mz']
        ms2_rt = ident_info_se['MS2_scan_time']
        dda = ident_info_se['DDA#']
        ms2_scan_id = ident_info_se['Scan#']
        ident_abbr = ident_info_se['Proposed_structures']
        ident_abbr = ident_abbr.replace('<', '&lt;')
        ident_abbr = ident_abbr.replace('>', '&gt;')
        score = ident_info_se['Overall_score']
        formula_ion = ident_info_se['Formula_ion']
        charge = ident_info_se['Charge']
        ident_idx = str(ident_idx)

        ident_info_df = pd.DataFrame()
        ident_info_df = ident_info_df.append(ident_info_se, ignore_index=True)

        with open(self.image_lst_page, 'a') as img_page:
            # convert info df to html table code
            if self.lipid_type in ['PA', 'PC', 'PE', 'PG', 'PI', 'PIP', 'PS']:
                plot_df_cols = ['Proposed_structures', 'Overall_score', 'i_sn1', 'i_sn2', 'i_[M-H]-sn1', 'i_[M-H]-sn2',
                                'i_[M-H]-sn1-H2O', 'i_[M-H]-sn2-H2O', 'SN_ratio']
            elif self.lipid_type in ['TG', 'TAG', 'DG', 'DAG', 'MG', 'MAG']:
                plot_df_cols = ['Proposed_structures', 'Overall_score', 'i_sn1', 'i_sn2', 'i_sn3',
                                'i_M-sn1', 'i_M-sn2', 'i_M-sn3',
                                'i_M-(sn1+sn2)', 'i_M-(sn1+sn3)', 'i_M-(sn2+sn3)', 'SN_ratio']
            else:
                plot_df_cols = ['Proposed_structures', 'Overall_score', 'i_sn1', 'i_sn2', 'i_[M-H]-sn1', 'i_[M-H]-sn2',
                                'i_[M-H]-sn1-H2O', 'i_[M-H]-sn2-H2O', 'SN_ratio']

            ident_col = ident_info_df.columns.tolist()

            for _col in plot_df_cols:
                if _col not in ident_col:
                    ident_info_df.loc[:, _col] = ''
            table_buf_code = ident_info_df.to_html(columns=plot_df_cols, float_format='%.1f', border=0, index=False)
            table_buf_code = table_buf_code.replace('NaN', '')
            img_title_str = ('{mz}_RT{rt:.3}_DDArank{dda}_Scan{scan}_{ident}_{f}_{chg}_score{score}'
                             .format(mz='%.4f' % ms1_pr_mz, rt=ms2_rt, dda=dda, scan=ms2_scan_id,
                                     ident=ident_abbr, score=score, f=formula_ion, chg=charge))
            img_info_lst = ['<a name="', ident_idx, '"><h3>', '<a href="', img_path, '" target="blank">', img_title_str,
                            '</a></h3></a>', '<a href="', img_path, '" target="blank">',
                            '<img src="', img_path, '" height="800" /></a>', table_buf_code, '\n<hr>\n']
            img_page.write(''.join(img_info_lst))

        with open(self.idx_lst_page, 'a') as idx_page:

            idx_str = ('''
                        <tr>\n
                        <td><a href ="LPPtiger_Results_Figures_list.html#{id}" target ="results_frame">{id}</td>\n
                        <td><a href ="LPPtiger_Results_Figures_list.html#{id}" target ="results_frame">{mz}</td>\n
                        <td><a href ="LPPtiger_Results_Figures_list.html#{id}" target ="results_frame">{rt}</td>\n
                        <td><a href ="LPPtiger_Results_Figures_list.html#{id}" target ="results_frame">{ident}</td>\n
                        <td><a href ="LPPtiger_Results_Figures_list.html#{id}" target ="results_frame">{score}</td>\n
                        </tr>\n
                        '''.format(id=ident_idx, mz='%.4f' % ms1_pr_mz, rt='%.1f' % ms2_rt,
                                   ident=ident_abbr, score=score))
            idx_page.write(idx_str)

        print('==> info added to report html -->')

    def add_all_info(self, ident_info_df):

        with open(self.image_lst_page, 'a') as img_page:
            with open(self.idx_lst_page, 'a') as idx_page:

                for _idx, _row_se in ident_info_df.iterrows():
                    img_path = _row_se['img_name']
                    ms1_pr_mz = _row_se['MS1_obs_mz']
                    ms2_rt = _row_se['MS2_scan_time']
                    dda = _row_se['DDA#']
                    ms2_scan_id = _row_se['Scan#']
                    ident_abbr = _row_se['Proposed_structures']
                    ident_abbr = ident_abbr.replace('<', '&lt;')
                    ident_abbr = ident_abbr.replace('>', '&gt;')
                    score = _row_se['Overall_score']
                    formula_ion = _row_se['Formula_ion']
                    charge = _row_se['Charge']
                    # ident_idx = str(_idx)

                    ident_info_df = pd.DataFrame()
                    ident_info_df = ident_info_df.append(_row_se, ignore_index=True)

                    # convert info df to html table code
                    if self.lipid_type in ['PA', 'PC', 'PE', 'PG', 'PI', 'PIP', 'PS']:
                        plot_df_cols = ['Proposed_structures', 'Overall_score', 'i_sn1', 'i_sn2', 'i_[M-H]-sn1',
                                        'i_[M-H]-sn2',
                                        'i_[M-H]-sn1-H2O', 'i_[M-H]-sn2-H2O', 'SN_ratio']
                    elif self.lipid_type in ['TG', 'TAG', 'DG', 'DAG', 'MG', 'MAG']:
                        plot_df_cols = ['Proposed_structures', 'Overall_score', 'i_sn1', 'i_sn2', 'i_sn3',
                                        'i_M-sn1', 'i_M-sn2', 'i_M-sn3',
                                        'i_M-(sn1+sn2)', 'i_M-(sn1+sn3)', 'i_M-(sn2+sn3)', 'SN_ratio']
                    else:
                        plot_df_cols = ['Proposed_structures', 'Overall_score', 'i_sn1', 'i_sn2', 'i_[M-H]-sn1',
                                        'i_[M-H]-sn2',
                                        'i_[M-H]-sn1-H2O', 'i_[M-H]-sn2-H2O', 'SN_ratio']

                    ident_info_df = pd.DataFrame(ident_info_df, columns=plot_df_cols)
                    ident_col = ident_info_df.columns.tolist()

                    for _col in plot_df_cols:
                        if _col not in ident_col:
                            ident_info_df.loc[:, _col] = 0
                    try:
                        table_buf_code = ident_info_df.to_html(float_format='%.1f', border=0, index=False)
                    except TypeError:
                        table_buf_code = ident_info_df.to_html(index=False)
                    table_buf_code = table_buf_code.replace('NaN', '')
                    img_title_str = ('{mz}_RT{rt:.3}_DDArank{dda}_Scan{scan}_{ident}_{f}_{chg}_score{score}'
                                     .format(mz='%.4f' % ms1_pr_mz, rt=ms2_rt, dda=dda, scan=ms2_scan_id,
                                             ident=ident_abbr, score=score, f=formula_ion, chg=charge))
                    img_info_lst = ['<a name="', '%i' % _idx, '"><h3>', '<a href="', img_path, '" target="blank">',
                                    img_title_str,
                                    '</a></h3></a>', '<a href="', img_path, '" target="blank">',
                                    '<img src="', img_path, '" height="800" /></a>', table_buf_code, '\n<hr>\n']
                    img_page.write(''.join(img_info_lst))

                    idx_str = ('''
                            <tr>\n
                            <td>
                            <a href ="LPPtiger_Results_Figures_list.html#{id}" target ="results_frame">{id}
                            </td>\n<td>
                            <a href ="LPPtiger_Results_Figures_list.html#{id}" target ="results_frame">{mz}
                            </td>\n<td>
                            <a href ="LPPtiger_Results_Figures_list.html#{id}" target ="results_frame">{rt}
                            </td>\n<td>
                            <a href ="LPPtiger_Results_Figures_list.html#{id}" target ="results_frame">{ident}
                            </td>\n<td>
                            <a href ="LPPtiger_Results_Figures_list.html#{id}" target ="results_frame">{score}
                            </td>\n</tr>\n
                            '''.format(id='%i' % _idx, mz='%.4f' % ms1_pr_mz, rt='%.1f' % ms2_rt,
                                       ident=ident_abbr, score=score))
                    idx_page.write(idx_str)

            print('==> info added to report html -->')

    def close_page(self):
        with open(self.main_page, 'a') as _m_page:
            _m_page.write('\n</body></html>\n')

        with open(self.image_lst_page, 'a') as _img_page:
            _img_page.write('\n</body></html>\n')

        with open(self.idx_lst_page, 'a') as _idx_page:
            _idx_page.write('\n</tbody>\n</table>\n</body></html>\n')
