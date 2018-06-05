import os, sys, re, json, string, random
from annotate_bed import MainApp
from flask import Flask, render_template, request, redirect, flash, send_from_directory, session, jsonify, copy_current_request_context
from flask_wtf import FlaskForm
from flask_wtf.file import FileAllowed, FileField, FileRequired
from werkzeug.utils import secure_filename
from wtforms import StringField, SubmitField, FileField, RadioField, SelectMultipleField
from wtforms.widgets import ListWidget, CheckboxInput
from wtforms.validators import DataRequired, Length, Required
from flask import url_for
from flask_bootstrap import Bootstrap
from multiprocessing.pool import ThreadPool
import re

from pymongo import MongoClient

app = Flask(__name__)
bootstrap = Bootstrap(app)
pool = ThreadPool(processes=1)
finished = False

UPLOAD_FOLDER = 'static/uploads'
ALLOWED_EXTENSIONS = {'txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif'}
ANNOT_CHOICES = [('all_beds','All'), ('circRNA','circRNA'), ('genes', 'Genes'),('coding_gene','Coding genes'),
                 ('longNC','Long non-coding'), ('mirna','miRNA'), ('mirbase','mirBase'),
                 ('noncoding_gene','Non-coding genes'),('pseudogene','Pseudogenes')]

GENELISTS = [('all_genes','All'), ('ASD', 'ASD'), ('ID_a', 'ID a'), ('ID_b','ID b'), ('dosage_sensitive', 'Dosage sensitive'),
             ('epilessia', 'Epilepsy'), ('malformazioni', 'Malformations'), ('mendeliome', 'Mendeliome'),
             ('onologhi', 'Onologs'), ('pubmed_autism_09-02-2018', 'PubMed Autism'),
             ('pubmed_brain_malformations_09-02-2018', 'PubMed Brain Malformations'),
             ('pubmed_epilepsy_or_seizures_09-02-2018', 'PubMed Epilepsy'),
             ('pubmed_intellectual_disability_09-02-2018', 'PubMed Intellectual Disability')]

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
def id_generator(size=8, chars=string.ascii_lowercase + string.digits):
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
    line_input = StringField(u'CNV string ', validators=[Length(max=50)])
    upload = FileField('Input file', validators=[
        FileAllowed(['txt', 'csv', 'cnv'], 'Text only!')
    ])
    window = StringField(u'Window (bp):', validators=[Length(max=15)], default='1000000')
    annot = MultiCheckboxField('Label', choices=ANNOT_CHOICES)
    genes = MultiCheckboxField('Label', choices=GENELISTS)
    submit = SubmitField("Submit")



# class FileForm(FlaskForm):
#     filein = FileField("Input file", validators=[DataRequired()])
#     submit = SubmitField("Submit")

@app.route('/', methods=['GET', 'POST'])
def index():
    form = MainForm()

    if form.validate_on_submit():
        

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
            session['working_filename'] = os.path.join(session['task_id'], session['task_id']+'.csv')
            session['cnv_line'] = request.form['line_input']
        return redirect(url_for('working'))

    return render_template('index.html', form=form)


@app.route('/working.html', methods=['GET', 'POST'])
def working():
    global async_result
    global finished
    finished = False

    args = Object
    session['file_out'] = os.path.join(app.config['UPLOAD_FOLDER'], "{}.xlsx".format(
        os.path.splitext(session['working_filename'])[0]))
    if session['choice'] == 'file':
        session['download_name'] = os.path.splitext(session['filename'])[0] + '_CANVAS.xlsx'
        setattrs(args, cnv_line=None)
        setattrs(args, cnv_file=os.path.join(app.config['UPLOAD_FOLDER'], session['working_filename']))
    elif session['choice'] == 'line':
        session['download_name'] = 'CANVAS_results.xlsx'
        setattrs(args, cnv_line=session['cnv_line'])
        setattrs(args, cnv_file=None)

    session['ann_choices'] = ['gene', 'mirna', 'pseudogene', 'longNC', 'circRNA'] #todo provvisori ovviamente
    setattrs(args, all=True, circRNA=False,
             coding_gene=False, distance=1000000, gene=False, longNC=False, mirna=False, mirbase=False,
             noncoding_gene=False, out=session['file_out'],
             pseudogene=False)

       
    async_result = pool.apply_async(worker, (args,))
    
    return render_template('working.html')
    # return render_template('results.html', file_out=file_out, download_name=download_name, result_db=result_db)


def worker(args):
    global finished
    result_db = MainApp(args).process()
    finished = True
    
@app.route('/status')
def thread_status():
    global async_result
    """ Return the status of the worker thread """
    return jsonify(dict(status=('finished' if finished else 'running')))

@app.route('/results.html', methods=['GET', 'POST'])
def results():
    return render_template('results.html', json_out=re.sub('.xlsx', '.json', session['file_out']),
                           file_out=session['file_out'], download_name=session['download_name'],
                           choices=session['ann_choices'])

if __name__ == '__main__':
    app.run()
