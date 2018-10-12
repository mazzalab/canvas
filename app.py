import os, sys, re, json, string, random
from annotate_bed import MainApp
from overlapper import OverlapApp
from flask import Flask, render_template, request, redirect, flash, send_from_directory, session, jsonify, copy_current_request_context
from flask_wtf import FlaskForm
from flask_wtf.file import FileAllowed, FileField, FileRequired
from werkzeug.utils import secure_filename
from wtforms import StringField, SubmitField, FileField, RadioField, SelectMultipleField
from wtforms.widgets import ListWidget, CheckboxInput
from wtforms.validators import DataRequired, Length, Required, regexp, data_required
from flask import url_for
from flask_bootstrap import Bootstrap
from multiprocessing.pool import ThreadPool
import re
import glob
from urllib.request import urlopen, urlretrieve
import json

from pymongo import MongoClient

app = Flask(__name__)
bootstrap = Bootstrap(app)
pool = ThreadPool(processes=1)
finished = False
UPLOAD_FOLDER = 'static/uploads'
ALLOWED_EXTENSIONS = {'txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif'}
ANNOT_CHOICES = [('all_beds','All'), ('coding_gene','Coding genes'),
                 ('noncoding_gene','Non-coding genes'), ('gene', 'Gene lists'), ('longNC','Long non-coding'),
                 ('mirna','MicroRNAs (miRNAs)'),
                 ('pseudogene','Pseudogenes'), ('circRNA','CircularRNAs (circRNAs)'), ('enhancer','Enhancers'),
                 ('ucr', 'Ultra Conserved Regions (UCRs)'), ('har', 'Human Accelerated Regions (HARs)')]


EXPLANATIONS = {"gene": "All RefSeq genes reported in UCSC genome browser.",
                "coding_gene": "RefSeq protein-coding genes (labelled as NM_*) reported in UCSC genome browser.",
                "noncoding_gene": "RefSeq non-protein-coding genes (labelled as NR_*) reported in UCSC genome browser.",
                "longNC": "Long non-coding RNAs reported in LNCipedia.",
                "mirna": "microRNAs, reported in miRBase, and microRNA targets, reported in DIANA-TarBase and "
                        "in TargetScan.",
                "circRNA": "Circular RNAs reported in circBase.", "pseudogene": "Pseudogenes reported in psiDR.",
                "ucr": "Ultraconserved elements (UCRs) reported in UCbase.",
                "har": "Human Accelerated Regions (HARs).",
                "enhancer": "Enhancers reported in Human Enhancer Disease Database (HEDD).",
                "ID_genelist": "Genes reported to be associated with intellectual disability (isolated and "
                               "associated disorders) in Vissers et al., 2016",
                "dosage_sensitive_genelist": "Genes reported to be dosage sensitive according to ClinGen "
                                             "Dosage Sensitivity Map.",
                "mendeliome_genelist": "Genes associated to Mendelian diseases reported in TruSight One "
                                       "Gene List (2013)",
                "ohnologs_genelist": "Ohnolog genes reported in Makino and McLysaght, 2010",
                "imprinted_genelist": "Mammalian imprinted genes taken from Catalogue of Parent of Origin "
                                      "Effects (2016)."
                }

EXPLANATIONS_TAD = {}

NICE_NAMES = {"gene": "Gene lists",
                "coding_gene": "Coding genes",
                "noncoding_gene": "Non-coding genes",
                "longNC": "Long NC",
                "mirna": "MiRNAs",
                "mirbase": "mirRBase",
                "circRNA": "CircRNAs",
                "pseudogene": "Pseudogenes",
                "ucr": "UCRs",
                "har": "HARs",
                "enhancer": "Enhancers",
                "ID_genelist": "Intellectual Disability",
                "dosage_sensitive_genelist": "Dosage sensitive",
                "mendeliome_genelist": "Mendeliome",
                "ohnologs_genelist":"Ohnolog",
                "imprinted_genelist": "Imprinted"
                }

NICE_NAMES_TAD = {}

GENELISTS = [('all_genelists','All'), ('ID', 'Intellectual Disability'), ('dosage_sensitive', 'Dosage sensitive'),
             ('mendeliome', 'Mendeliome panel'),
             ('ohnologs', 'Ohnologs'), ('imprinted', 'Imprinted')]

SOURCES = {
    "gene": "<a href='https://genome.ucsc.edu/cgi-bin/hgTables' target='_blank'>UCSC</a>",
    "coding_gene": "<a href='https://genome.ucsc.edu/cgi-bin/hgTables' target='_blank'>UCSC</a>",
    "noncoding_gene": "<a href='https://genome.ucsc.edu/cgi-bin/hgTables' target='_blank'>UCSC</a>",
    "longNC": "<a href='https://lncipedia.org/' target='_blank'>LNCipedia</a>",
    "mirna": "<a href='ftp://mirbase.org/pub/mirbase/20/genomes/' target='_blank'>miRBase</a>",
    "mirbase": "<a href='ftp://mirbase.org/pub/mirbase/20/genomes/' target='_blank'>miRBase</a>",
    "circRNA": "<a href='http://circbase.org/cgi-bin/downloads.cgi' target='_blank'>circBase</a>",
    "pseudogene": "<a href='http://www.pseudogenes.org/psidr/' target='_blank'>psiDR</a>",
    "ucr": "<a href='http://ucbase.unimore.it/' target='_blank'>UCbase</a>",
    "har": "<a href='https://www.ncbi.nlm.nih.gov/pubmed/27667684' target='_blank'>Doan et al., 2016</a>",
    "enhancer": "<a href='http://zdzlab.einstein.yu.edu/1/hedd/download.php' target='_blank'>HEDD</a>",
    "ID_genelist": "<a href='https://www.ncbi.nlm.nih.gov/pubmed/26503795' target='_blank'>Vissers et al., 2016</a>",
    "dosage_sensitive_genelist": "<a href='https://www.ncbi.nlm.nih.gov/projects/dbvar/clingen/help.shtml' target='_blank'>NCBI</a>",
    "mendeliome_genelist": "<a href='https://support.illumina.com/downloads/trusight_one_sequencing_panel_product_file.html' target='_blank'>Illumina.com</a>",
    "ohnologs_genelist": "<a href='http://www.pnas.org/content/107/20/9270' target='_blank'>Makino and McLysaght, 2010</a>",
    "imprinted_genelist": "<a href='http://igc.otago.ac.nz/1601summarytable.pdf' target='_blank'>Catalogue of Parent of Origin Effects</a>",
}

#it's value, label
TISSUE_CHOICES = [('all', 'All'), ('Liver_STL011_Leung2015', 'Liver_STL011_Leung2015'),
               ('Lung_Donor_LG1', 'Lung_Donor_LG1'),
               ('LNCaP_rep1', 'LNCaP_rep1'), ('Pancreas_Donor_PA2', 'Pancreas_Donor_PA2'),
               ('SKNMC_rep1', 'SKNMC_rep1'), ('G401_rep1', 'G401_rep1'),
               ('H1_ESC_Dixon2015', 'H1_ESC_Dixon2015'), ('GM12878_Lieberman', 'GM12878_Lieberman'),
               ('H1_NPC_Dixon2015', 'Long description'), ('K562_Lieberman', 'Long description'),
               ('SKMEL5_rep1', 'Long description'), ('SKNDZ_rep1', 'Long description'),
               ('AdrenalGland_Donor_AD2', 'Long description'), ('Caki2_rep1', 'Long description'),
               ('T470_rep1', 'Long description'), ('SJCRH30_rep1', 'Long description'),
               ('HMEC_Lieberman', 'Long description'), ('NCIH460_rep1', 'Long description'),
               ('PANC1_rep1', 'Long description'), ('H1_TRO_Dixon2015', 'Long description'),
               ('HUVEC_Lieberman', 'Long description'), ('Aorta_STL002_Leung2015', 'Long description'),
               ('Cortex_DLPFC_Donor_CO', 'Long description'),
               ('VentricleRight_Donor_RV3', 'Long description'),
               ('VentricleLeft_STL003_Leung2015', 'Long description'),
               ('Spleen_Donor_PX1', 'Long description'), ('H1_MSC_Dixon2015', 'Long description'),
               ('A549_rep1', 'Long description'), ('Bladder_Donor_BL1', 'Long description'),
               ('H1_MES_Dixon2015', 'Long description'), ('Bowel_Small_Donor_SB2', 'Long description'),
               ('IMR90_Lieberman', 'Long description'), ('NHEK_Lieberman', 'Long description'),
               ('RPMI7951_rep1', 'Long description'), ('KBM7_Lieberman', 'Long description'),
               ('Muscle_Psoas_Donor_PO1', 'Long description'), ('Thymus_STL001_Leung2015', 'Long description')]

TISSUE_DESCRIPTIONS = {'All': 'All availabe tissues', 'Liver_STL011_Leung2015': 'Loooooooooooooooooooooooong description',
                       'Lung_Donor_LG1': 'Long description',
                       'LNCaP_rep1': 'Long description', 'Pancreas_Donor_PA2': 'Long description',
                       'SKNMC_rep1': 'Long description', 'G401_rep1': 'Long description',
                       'H1_ESC_Dixon2015': 'Long description', 'GM12878_Lieberman': 'Long description',
                       'H1_NPC_Dixon2015': 'Long description', 'K562_Lieberman': 'Long description',
                       'SKMEL5_rep1': 'Long description', 'SKNDZ_rep1': 'Long description',
                       'AdrenalGland_Donor_AD2': 'Long description', 'Caki2_rep1': 'Long description',
                       'T470_rep1': 'Long description', 'SJCRH30_rep1': 'Long description',
                       'HMEC_Lieberman': 'Long description', 'NCIH460_rep1': 'Long description',
                       'PANC1_rep1': 'Long description', 'H1_TRO_Dixon2015': 'Long description',
                       'HUVEC_Lieberman': 'Long description', 'Aorta_STL002_Leung2015': 'Long description',
                       'Cortex_DLPFC_Donor_CO': 'Long description',
                       'VentricleRight_Donor_RV3': 'Long description',
                       'VentricleLeft_STL003_Leung2015': 'Long description',
                       'Spleen_Donor_PX1': 'Long description', 'H1_MSC_Dixon2015': 'Long description',
                       'A549_rep1': 'Long description', 'Bladder_Donor_BL1': 'Long description',
                       'H1_MES_Dixon2015': 'Long description', 'Bowel_Small_Donor_SB2': 'Long description',
                       'IMR90_Lieberman': 'Long description', 'NHEK_Lieberman': 'Long description',
                       'RPMI7951_rep1': 'Long description', 'KBM7_Lieberman': 'Long description',
                       'Muscle_Psoas_Donor_PO1': 'Long description', 'Thymus_STL001_Leung2015': 'Long description'}

app.config['SECRET_KEY'] = 'AGATTAcanvas2018'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_PATH'] = 5000



#
#     return '''
#     <!doctype html>
#     <title>Upload new File</title>
#     <h1>Upload new File</h1>
#     <form method=post enctype=multipart/form-data>
#       <p><input type=file name=file>
#          <input type=submit value=Upload>
#     </form>
#     '''
#
# @app.route('/uploads/<filename>')
# def uploaded_file(filename):
#     return send_from_directory(app.config['UPLOAD_FOLDER'],
#
#          filename)

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def id_generator(size=12, chars=string.ascii_lowercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


class Object(object):
    pass


def setattrs(_self, **kwargs):
    for k,v in kwargs.items():
        setattr(_self, k, v)
        
class MultiCheckboxField(SelectMultipleField):
    widget = ListWidget(prefix_label=False)
    option_widget = CheckboxInput()
    
class MainForm(FlaskForm):
    line_input = StringField(u'Genomic region ', validators=[Length(max=50)])
    upload = FileField('Input file', validators=[
        FileAllowed(['txt', 'csv', 'cnv'], 'Text only!')
    ])
    
    overlap_upload = FileField('Input file', validators=[
        FileAllowed(['txt', 'csv', 'cnv'], 'Text only!')
    ])
    window = StringField(u'Window (bp):', validators=[Length(max=15)], default='1000000')
    window_tad = StringField(u'Window (bp):', validators=[Length(max=15)], default='1000000')
    padding = StringField(u'Padding (bp):', validators=[Length(max=15)], default='0')
    min_ovl_rec = StringField(u'Min. Overlap (%):', validators=[Length(max=5)], default='50')
    min_ovl_span = StringField(u'Min. Overlap (%):', validators=[Length(max=5)], default='50')
    max_span = StringField(u'Max. Span (bp):', validators=[Length(max=15)], default='100000')
    annot = MultiCheckboxField('annot', choices=ANNOT_CHOICES)
    genes = MultiCheckboxField('genes', choices=GENELISTS)
    # tissues = SelectMultipleField('Tissues', choices=TISSUE_CHOICES)
    # print(tissues)
    submit = SubmitField("Submit")



# class FileForm(FlaskForm):
#     filein = FileField("Input file", validators=[DataRequired()])
#     submit = SubmitField("Submit")

@app.route('/', methods=['GET', 'POST'])
def index():
    session['ann_choices'] = []
    session['genes_choices'] = []
    session['choice'] = ''
    session['cnv_line'] = None
    session['filename'] = ''
    session['working_filename'] = ''
    session['task_id'] = ''
    session['window'] = []
    session['file_out'] = ''
    session['ref'] = ''
    session['padding'] = ''
    session['min_ovl_rec'] = ''
    session['min_ovl_span'] = ''
    session['max_span'] = ''
    session['interset_choice'] = ''
    session['intraset_choice'] = False
    session['ovl_file'] = ''
    session['working_ovl_filename'] = ''
    session['overlap_upload'] = ''
    session['combine_mode'] = ''
    session['tissue_choices'] = []
    session['window_tad'] = ''
    
    
    # session['overlap_fileout_xlsx'] = ''
    # session['overlap_fileout_csv'] =
    
    print("SESSIONE")
    print(session['cnv_line'])
    print(session['task_id'])
    form = MainForm()

    if form.validate_on_submit():
        print("REQUEST FORM:")
        print(request.form)
        session['task_id'] = id_generator()
        os.mkdir(os.path.join(app.config['UPLOAD_FOLDER'], session['task_id']))
        if 'radio2' in request.form:
            session['choice'] = 'file'
            session['cnv_line'] = None
            f = form.upload.data
            session['filename'] = secure_filename(f.filename)
            session['working_filename'] = os.path.join(session['task_id'], session['task_id']+'.csv')
            f.save(os.path.join(
                app.config['UPLOAD_FOLDER'], session['working_filename']
            ))
        elif 'radio1' in request.form:
            session['choice'] = 'line'
            session['filename'] = 'line input'
            session['working_filename'] = os.path.join(session['task_id'], session['task_id']+'.csv')
            session['cnv_line'] = request.form['line_input']

        
        if 'radio-hg19' in request.form:
            session['ref'] = 'hg19'
        elif 'radio-hg18' in request.form:
            session['ref'] = 'hg18'
        elif 'radio-hg38' in request.form:
            session['ref'] = 'hg38'
        print("REF:", session['ref'])

        if 'annotradio' in request.form:
            for elem in request.form.getlist('annot'):
                    session['ann_choices'].append(elem)
            if 'mirna' in request.form.getlist('annot'):
                session['ann_choices'].append('mirbase')
            
            for elem in request.form.getlist('genes'):
                    if elem != "all_genelists":
                        session['genes_choices'].append(elem+'_genelist')
                
            if 'window' in request.form:
                session['window'] = request.form['window']
            else:
                session['window'] = 1000000
                
            print("ANNOT:", session['ann_choices'])
            print("GENES:", session['genes_choices'])
            print("DISTANCE", session['window'])
            return redirect(url_for('working'))
        
        elif 'ovlradio' in request.form:
            if 'line_input' in request.form:
                with open(os.path.join(app.config['UPLOAD_FOLDER'], session['working_filename']), 'w') as f:
                    m = re.match(r'(?P<chr>chr[\dXYM]+):(?P<start>\d+)-(?P<end>\d+)', session['cnv_line'])
                    
                    f.write("CHR\tSTART\tEND\n{0}\t{1}\t{2}".format(m.group('chr'), m.group('start'),
                                                                    m.group('end')))
                
            #Intraset or interset
            if 'radio_interset' in request.form:
                session['interset_choice'] = request.form['overlapselect']
                if request.form['overlapselect'] == 'FILE':
                    f = form.overlap_upload.data
                    session['overlap_upload'] = secure_filename(f.filename)
                    session['working_ovl_filename'] = os.path.join(session['task_id'],
                                                               session['task_id'] + '_ovl.csv')
                    f.save(os.path.join(
                        app.config['UPLOAD_FOLDER'], session['working_ovl_filename']
                    ))
                    session['ovl_file'] = os.path.join(app.config['UPLOAD_FOLDER'],
                                                        session['working_ovl_filename'])
                    
                elif request.form['overlapselect'] == 'DGV_overlap':
                    print("DGV CHOSEN!!!!!!!!!!!!!")
                    session['ovl_file'] = 'DGV_overlap'
                    
            elif 'radio_intraset' in request.form:
                session['intraset_choice'] = True
                session['ovl_file'] = os.path.join(app.config['UPLOAD_FOLDER'], session['working_filename'])

            # Reciprocal or spanning
            if 'radio_reciprocal' in request.form:
                session['ovl_mode'] = 'reciprocal'
                session['padding'] = int(request.form['padding'])
                session['min_ovl_rec'] = int(request.form['min_ovl_rec'])
                
            elif 'radio_spanning' in request.form:
                session['ovl_mode'] = 'spanning'
                session['min_ovl_span'] = int(request.form['min_ovl_span'])
                session['max_span'] = int(request.form['max_span'])
            
            print("Sessione", session.__dict__)
            return redirect(url_for('working_ovl'))
    
        elif 'tadradio' in request.form:
            for elem in request.form.getlist('tad-tissue'):
                    if elem == 'All':
                        session['tissue_choices'] = list(TISSUE_DESCRIPTIONS.keys())
                        session['tissue_choices'].remove('All')
                        break
                    else:
                        session['tissue_choices'].append(elem)
            
            session['tissue_choices'] = ','.join(session['tissue_choices'])
            
            if 'window_tad' in request.form:
                session['window_tad'] = request.form['window_tad']
            else:
                session['window_tad'] = 1000000
                
            print(session['tissue_choices'])
            return redirect(url_for('working_tad'))

        

    return render_template('index.html', form=form, tissue_descr=TISSUE_DESCRIPTIONS)


@app.route('/working.html', methods=['GET', 'POST'])
def working():
    global async_result
    global finished
    finished = False

    session['file_out'] = os.path.join(app.config['UPLOAD_FOLDER'], "{}.xlsx".format(
        os.path.splitext(session['working_filename'])[0]))
    
    args = Object()
    print("ARGOMENTI BEFORE")
    print(args.__dict__)
    setattrs(args, cnv_file=None, cnv_line=None, all_beds=False, circRNA=False,
             coding_gene=False, gene=False, longNC=False, mirna=False, mirbase=False,
             noncoding_gene=False, pseudogene=False, ucr=False, har=False, enhancer=False,
             all_genelists=False,
             ID_genelist=False, dosage_sensitive_genelist=False, imprinted_genelist=False,
             mendeliome_genelist=False,
             ohnologs_genelist=False, distance=int(session['window']), reference='hg19',
             out=session['file_out'])
    
    if session['choice'] == 'file':
        session['download_name'] = os.path.splitext(session['filename'])[0] + '_INCAS.xlsx'
        setattrs(args, cnv_line=None)
        setattrs(args, cnv_file=os.path.join(app.config['UPLOAD_FOLDER'], session['working_filename']))
    elif session['choice'] == 'line':
        session['download_name'] = 'INCAS_results.xlsx'
        setattrs(args, cnv_line=session['cnv_line'])
        setattrs(args, cnv_file=None)
    
    if session['ref'] != 'hg19':
        setattr(args, 'reference', session['ref'])
        
    for elem in session['ann_choices']:
        setattr(args, elem, True)
    
    for elem in session['genes_choices']:
        setattr(args, elem, True)
  

    print("GLI ARGOMENTI")
    print(args.__dict__)
    async_result = pool.apply_async(worker, (args,))
    
    # Comic load!
    comicnum = random.randint(1,2017)
    url = 'https://xkcd.com/' + str(comicnum) + '/info.0.json'
    u = urlopen(url)
    page_html = u.read()
    u.close()

    json_data = json.loads(page_html)
    session['comic'] = json_data['img']

    print("FINISHED the working")
    return render_template('working.html', file=session['filename'], nice_names=NICE_NAMES, comic=session['comic'])
    # return render_template('results.html', file_out=file_out, download_name=download_name, result_db=result_db)


@app.route('/working_ovl.html', methods=['GET', 'POST'])
def working_ovl():
    global async_result
    global finished
    finished = False
    
    session['file_out'] = os.path.join(app.config['UPLOAD_FOLDER'], "{}".format(
        os.path.splitext(session['working_filename'])[0]+'_matrix.xlsx'))
    
    args = Object()
    print("ARGOMENTI BEFORE")
    print(args.__dict__)
    setattrs(args,  combine_mode=None,
             input_1=None, input_2=None, mode=None, min_overlap=None, output_prefix=None,
             padding=None, span=10000)

    setattrs(args, input_1=os.path.join(app.config['UPLOAD_FOLDER'], session['working_filename']))
    setattrs(args, input_2=session['ovl_file'])
    setattrs(args, output_prefix=os.path.join(app.config['UPLOAD_FOLDER'], session['task_id'], session['task_id']))
    
    if session['choice'] == 'file':
        session['download_name'] = os.path.splitext(session['filename'])[0] + '_INCAS_overlap_matrix.xlsx'

    elif session['choice'] == 'line':
        session['download_name'] = 'INCAS_overlap_matrix.xlsx'

    # if session['ref'] != 'hg19':
    #     setattr(args, 'reference', session['ref'])
    
    if session['intraset_choice'] is True:
        setattrs(args, combine_mode='combination')
    elif session['interset_choice'] != '':
        setattrs(args, combine_mode='product')

    if session['ovl_mode'] == 'reciprocal':
        setattrs(args, mode='reciprocal')
        setattrs(args, min_overlap=session['min_ovl_rec'])
        setattrs(args, padding=session['padding'])
        
    elif session['ovl_mode'] == 'spanning':
        setattrs(args, mode='spanning')
        setattrs(args, span=session['max_span'])
        setattrs(args, min_overlap=session['min_ovl_span'])
        setattrs(args, padding=0)


    print("GLI ARGOMENTI")
    print(args.__dict__)

    async_result = pool.apply_async(worker_ovl, (args,))
    
    # Comic load!
    comicnum = random.randint(1, 2017)
    url = 'https://xkcd.com/' + str(comicnum) + '/info.0.json'
    u = urlopen(url)
    page_html = u.read()
    u.close()
    
    json_data = json.loads(page_html)
    session['comic'] = json_data['img']
    
    print("FINISHED the working_ovl")
    return render_template('working_ovl.html', file=session['filename'], nice_names=NICE_NAMES,
                           comic=session['comic'])


@app.route('/working_tad.html', methods=['GET', 'POST'])
def working_tad():
    global async_result
    global finished
    finished = False
    
    session['file_out'] = os.path.join(app.config['UPLOAD_FOLDER'], "{}.xlsx".format(
        os.path.splitext(session['working_filename'])[0]))
    
    args = Object()
    print("ARGOMENTI BEFORE")
    print(args.__dict__)
    setattrs(args, cnv_file=None, cnv_line=None, TAD=None, distance=int(session['window_tad']), reference='hg19',
             out=session['file_out'])
    
    if session['choice'] == 'file':
        session['download_name'] = os.path.splitext(session['filename'])[0] + '_INCAS.xlsx'
        setattrs(args, cnv_line=None)
        setattrs(args, cnv_file=os.path.join(app.config['UPLOAD_FOLDER'], session['working_filename']))
    elif session['choice'] == 'line':
        session['download_name'] = 'INCAS_results.xlsx'
        setattrs(args, cnv_line=session['cnv_line'])
        setattrs(args, cnv_file=None)
    
    if session['ref'] != 'hg19':
        setattr(args, 'reference', session['ref'])
    

    setattr(args, 'TAD', session['tissue_choices'])
    

    print("GLI ARGOMENTI")
    print(args.__dict__)
    async_result = pool.apply_async(worker_tad, (args,))
    
    # Comic load!
    comicnum = random.randint(1, 2017)
    url = 'https://xkcd.com/' + str(comicnum) + '/info.0.json'
    u = urlopen(url)
    page_html = u.read()
    u.close()
    
    json_data = json.loads(page_html)
    session['comic'] = json_data['img']
    
    print("FINISHED the working")
    return render_template('working_tad.html', file=session['filename'], nice_names=NICE_NAMES,
                           comic=session['comic'])

def worker(args):
    global finished
    success = MainApp(args).process()
    if success != 0:
        finished = -1
    else:
        finished = True
 
 
def worker_ovl(args):
    global finished
    success = OverlapApp(args).process()
    print("SUCCESS OVL", success)
    if success != 0:
        finished = -1
    else:
        finished = True

def worker_tad(args):
    global finished
    success = MainApp(args).process()
    print("SUCCESS TAD", success)
    if success != 0:
        finished = -1
    else:
        finished = True
        
@app.route('/status')
def thread_status():
    global async_result
    
    """ Return the status of the worker thread """
    # f_read = open(session['file_out'].replace('.xlsx','_log.txt'), "r").readlines()
    # # st_results = os.stat(session['file_out'].replace('.xlsx','_log.txt'))
    # # st_size = st_results[6]
    # # f_read.seek(st_size)
    progress = [os.path.basename(x).split('.')[0] for x in glob.glob(os.path.dirname(session['file_out'])+'/*.progress')]
    progress.sort(key=natural_keys)
    print("STATUS")
    print(progress)
    if finished == True:
        return jsonify(dict(status='finished'))
    elif finished == -1:
        return jsonify(dict(status='problem'))
    else:
        return jsonify(dict(status=progress))


@app.route('/status_ovl')
def thread_status_ovl():
    global async_result
    
    """ Return the status of the worker thread """
    # f_read = open(session['file_out'].replace('.xlsx','_log.txt'), "r").readlines()
    # # st_results = os.stat(session['file_out'].replace('.xlsx','_log.txt'))
    # # st_size = st_results[6]
    # # f_read.seek(st_size)
    # progress = [os.path.basename(x).split('.')[0] for x in
    #             glob.glob(os.path.dirname(session['file_out']) + '/*.progress')]
    # progress.sort(key=natural_keys)
    # print("STATUS")
    # print(progress)
    if finished == True:
        return jsonify(dict(status='finished'))
    elif finished == -1:
        return jsonify(dict(status='problem'))
    else:
        return jsonify(dict(status='Wait'))


@app.route('/status_tad')
def thread_status_tad():
    global async_result
    
    """ Return the status of the worker thread """
    # f_read = open(session['file_out'].replace('.xlsx','_log.txt'), "r").readlines()
    # # st_results = os.stat(session['file_out'].replace('.xlsx','_log.txt'))
    # # st_size = st_results[6]
    # # f_read.seek(st_size)
    # progress = [os.path.basename(x).split('.')[0] for x in
    #             glob.glob(os.path.dirname(session['file_out']) + '/*.progress')]
    # progress.sort(key=natural_keys)
    # print("STATUS")
    # print(progress)
    if finished == True:
        return jsonify(dict(status='finished'))
    elif finished == -1:
        return jsonify(dict(status='problem'))
    else:
        return jsonify(dict(status='Wait'))
    
@app.route('/results.html', methods=['GET', 'POST'])
def results():
    print("IN RESULTS")
    if 'all_beds' in session['ann_choices']:
        print(session)
        print(session['ann_choices'])
        session['ann_choices'].remove('all_beds')
        
    if 'all_genelists' in session['ann_choices']:
        session['ann_choices'].remove('all_genelists')
        
    return render_template('results.html', json_out=re.sub('.xlsx', '.json', session['file_out']),
                           file_out=session['file_out'],
                           text_file_out=re.sub('.xlsx', '.csv', session['file_out']),
                           download_name=session['download_name'],
                           text_download_name=re.sub('.xlsx', '.csv', session['download_name']),
                           choices=session['ann_choices']+session['genes_choices'],
                           genes_choices=session['genes_choices'],
                           distance=session['window'],
                           info=EXPLANATIONS, nice_names=NICE_NAMES, sources=SOURCES)


@app.route('/results_ovl.html', methods=['GET', 'POST'])
def results_ovl():
    print("IN RESULTS")
    
    return render_template('results_ovl.html',
                           file_out=session['file_out'],
                           text_file_out=re.sub('_matrix.xlsx', '_list.csv', session['file_out']),
                           download_name=session['download_name'],
                           text_download_name=re.sub('_matrix.xlsx', '_list.csv', session['download_name'])
                           )


@app.route('/results_tad.html', methods=['GET', 'POST'])
def results_tad():
    print("IN RESULTS")
    print(session['download_name'])
    print(session['file_out'])

    return render_template('results_tad.html',
                           file_out=session['file_out'], json_out=re.sub('.xlsx', '.json', session['file_out']),
                           text_file_out=re.sub('.xlsx', '.csv', session['file_out']),
                           choices=session['tissue_choices'].split(','),
                           download_name=session['download_name'],
                           text_download_name=re.sub('.xlsx', '.csv', session['download_name'],),
                           distance=session['window_tad'],
                           info=EXPLANATIONS_TAD
                           )

@app.route('/error.html', methods=['GET', 'POST'])
def problem():
    return render_template('error.html', file=session['filename'])

@app.route('/privacy-and-cookie-policy.html')
def privacy():
    return render_template('privacy-and-cookie-policy.html')

if __name__ == '__main__':
    app.run(debug=True)
