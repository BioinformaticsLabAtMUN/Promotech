import numpy as np, pandas as pd, time
from io import StringIO
from Bio import SeqIO
import joblib, traceback, os, sys
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QRadioButton, QStackedWidget, QPlainTextEdit, QLabel, QFileDialog, QProgressBar, QTableView
from PyQt5 import QtWidgets, uic, QtCore
from PyQt5.QtCore import QRunnable, pyqtSlot, pyqtSignal, QObject, QAbstractTableModel
Qt = QtCore.Qt

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append("{}/../sequences".format(dir_path))
sys.path.append("{}/../utils".format(dir_path))
from process_sequences import predictSequencesFromString
from utils import print_fn

class PandasModel(QAbstractTableModel):
    def __init__(self, data, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self._data = data
    def rowCount(self, parent=None):
        return len(self._data.values)
    def columnCount(self, parent=None):
        return self._data.columns.size
    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:
                return QtCore.QVariant(str(
                    self._data.values[index.row()][index.column()]))
        return QtCore.QVariant()
    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return str(self._data.columns[section])

            if orientation == Qt.Vertical:
                return str(self._data.index[section])

class ParrallelJobSignals(QObject):
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    result = pyqtSignal(object)
    progress = pyqtSignal(int)

class ParrallelJob(QRunnable):
    def __init__(self, fn, *args, **kwargs):
        super(ParrallelJob, self).__init__()
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = ParrallelJobSignals()
#    @pyqtSlot()
    def run(self):
        try:
            result = self.fn(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)
        finally:
            self.signals.finished.emit()

class Promotech_UI(QWidget):
    def __init__(
        self, ui_path="form.ui",
        init_function=None,
        preprocess_seqs_fn=None,
        predict_seqs_fn=None,
        predict_gen_fn=None
    ):
        super(Promotech_UI, self).__init__()
        self.ui_path=ui_path
        self.model_file_name = 'RF-HOT.model'
        self.model_type = "RF"
        self.model_mode = "HOT"
        self.model_types = ["RF", ""]
        self.load_ui()
        self.preprocess_seqs_fn = preprocess_seqs_fn
        self.predict_seqs_fn = predict_seqs_fn
        self.predict_gen_fn    = predict_gen_fn
        # LOAD MODEL IN PARALLEL THREAD
        # new_job         = ParrallelJob(fn=self.load_model)
        # self.threadpool = QtCore.QThreadPool()
        # self.threadpool.start(new_job)
        # new_job.signals.result.connect(self.prepare_model)

    def load_ui(self):
        self.window = uic.loadUi(self.ui_path, self)
        # LOAD WIDGETS
        #self.genome_radio_button           = self.window.findChild(QRadioButton   , "genome_radio_button")
        self.sequences_radio_button        = self.window.findChild(QRadioButton   , "sequences_radio_button")
        self.console_radio_button          = self.window.findChild(QRadioButton   , "console_radio_button")
        self.results_radio_button          = self.window.findChild(QRadioButton   , "results_radio_button")
        self.sequences_stacked_widget      = self.window.findChild(QStackedWidget , "sequences_stacked_widget")
        self.sequences_text_edit           = self.window.findChild(QPlainTextEdit , 'sequences_text_edit')
        self.message_label                 = self.window.findChild(QLabel         , 'message_label')
        self.sequences_file_location_label = self.window.findChild(QLabel         , 'sequences_file_location_label')
        self.sequences_predict_push_button = self.window.findChild(QPushButton    , 'sequences_predict_push_button')
        self.sequences_browse_push_button  = self.window.findChild(QPushButton    , 'sequences_browse_push_button')
        self.progress_bar                  = self.window.findChild(QProgressBar   , 'progress_bar')
        self.console_text_edit             = self.window.findChild(QPlainTextEdit , 'console_text_edit')
        self.results_table_view            = self.window.findChild(QTableView     , 'results_table_view')
        # SET SLOTS & EVENTS
        #self.genome_radio_button.toggled.connect(self.set_genome_prediction_view)
        self.sequences_radio_button.toggled.connect(self.set_sequence_prediction_view)
        self.results_radio_button.toggled.connect(self.set_results_view)
        self.console_radio_button.toggled.connect(self.set_console_view)
        self.sequences_predict_push_button.clicked.connect(self.predict_sequences)
        self.sequences_browse_push_button.clicked.connect(self.browse_sequences_file)

    def lock_widget(self, widget, is_enabled):
        if(widget is None):
            print("ERROR LOCKING THE MAIN SCREEN WHILE LOADING. WIDGET IS NULL.", widget, is_enabled)
            return
        widget.setEnabled(is_enabled)
    def load_model(self):
        self.lock_widget(self.sequences_stacked_widget, False)
        self.progress_bar.setRange(0, 0)
        self.progress_bar.setValue(0)
        model_path = os.path.join(dir_path, "../models/", self.model_file_name ) 
        self.update_message_label("LOADING MODEL...")
        model = joblib.load(model_path)
        self.update_message_label("MODEL LOADED SUCCESSFULLY...")
        return model
    def prepare_model(self, model):
        self.lock_widget(self.sequences_stacked_widget, True)
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.model = model
        print(model)
        self.update_message_label("MODEL LOADED SUCCESSFULLY.")
    def load_model_error(self, error_data ):
        self.model = model
        self.update_message_label("MODEL LOADED WITH ERRORS: "+str(error_data))
    def set_sequence_prediction_view(self, trigger):
        if(trigger):
            self.sequences_stacked_widget.setCurrentIndex(0)
            # self.update_message_label("SET TO 40-NT SEQUENCES PREDICTION MODE.")
            self.progress_bar.setValue(0)
    def set_genome_prediction_view(self, trigger):
        if(trigger):
            self.sequences_stacked_widget.setCurrentIndex(1)
            # self.update_message_label("SET TO GENOME PREDICTION MODE.")
            self.progress_bar.setValue(0)
    def set_console_view(self, trigger):
        if(trigger):
            self.sequences_stacked_widget.setCurrentIndex(2)
            # self.update_message_label("SET TO CONSOLE VIEW.")
    def set_results_view(self, trigger):
        if(trigger):
            self.sequences_stacked_widget.setCurrentIndex(3)
            # self.update_message_label("SET TO RESULTS VIEW.")
    def update_results_table( self, chroms, seqs, y_pred , prepro_seqs):
        df = pd.DataFrame([], columns = ['chroms', 'sequences', 'preprocessed_sequences', 'predictions', 'observation'])
        df["chroms"] = chroms
        df["sequences"] = seqs
        df["preprocessed_sequences"] = [ "".join(map(str, row)) for row in prepro_seqs.values ]
        df["predictions"] = y_pred
        df["observation"] = [ "PROMOTER (>= 0.5)" if pred >= 0.5 else "NON-PROMOTER (< 0.5)" for pred in y_pred ]

        self.sequences_stacked_widget.setCurrentIndex(3)
        self.results_radio_button.setChecked(True)
        self.results_model = PandasModel(data=df, parent=self)
        self.results_table_view.setModel(self.results_model)
        header = self.results_table_view.horizontalHeader()
        header_width_formats = [
            QtWidgets.QHeaderView.ResizeToContents, QtWidgets.QHeaderView.ResizeToContents,
            QtWidgets.QHeaderView.Interactive, QtWidgets.QHeaderView.ResizeToContents,
            QtWidgets.QHeaderView.ResizeToContents
        ]
        # header_width_formats = [ QtWidgets.QHeaderView.Interactive ] * len(df.columns)
        for i in range(len(df.columns)):
            header.setSectionResizeMode(i, header_width_formats[i] )
            
    def start_multi_thread_prediction(self):
        seqs_content      = self.sequences_text_edit.toPlainText()
        self.chroms, self.raw_seqs, self.y_pred, self.prepro_seqs = predictSequencesFromString(
            seqs_content, 
            print_fn=self.update_console, 
            out_dir="results", 
            threshold=0.5, 
            model_type="RF-HOT" , 
            log_file="./results/gui_prediction.log.txt", 
            start_time=time.time()
        )
        self.progress_bar.setRange(0, 100)
        # self.lock_widget(self.sequences_stacked_widget, True)
        self.update_results_table( self.chroms, self.raw_seqs, self.y_pred , self.prepro_seqs )
        # self.update_message_label( "PROMOTER PREDICTIONS SCORES GENERATED SUCCESSFULLY. X_INV THE RESULTS TAB FOR MORE INFORMATION.".format( len(self.raw_seqs) ) )
    def predict_sequences(self, print_fn=print):
        self.sequences_stacked_widget.setCurrentIndex(2)
        self.console_radio_button.setChecked(True)
        # self.lock_widget(self.sequences_stacked_widget, False)
        self.progress_bar.setRange(0, 0)
        new_job         = ParrallelJob(fn=self.start_multi_thread_prediction)
        self.threadpool = QtCore.QThreadPool()
        self.threadpool.start(new_job)
        
    def browse_sequences_file(self):
        self.sequences_file_path = QFileDialog.getOpenFileName(
            parent=self,
            directory="{}/../".format(dir_path),
            caption='Open FASTA file',
            filter="Image files (*.txt *.fasta *.fa)" )[0]
        print("BROWSE", self.sequences_file_path)
        self.sequences_file_location_label.setText(str(self.sequences_file_path))
        if(self.sequences_file_path == ''):
            return
        seqs_file         = open(self.sequences_file_path, "r")
        seqs_file_content = seqs_file.read()
        self.sequences_text_edit.setPlainText(seqs_file_content)

    def predict_genome(self):
        pass
    def update_message_label(self, msg):
        self.message_label.setText(msg)
    
    def update_console(self, msg, log_file ):
        print(msg)
        self.console_text_edit.insertPlainText( msg );


def demo_preprocess(seqs, progress_bar=None, print_fn=None):
    print("PREPROCESSING")
    preprocessed_seqs = parseSequences(seqs=seqs, sequence_length=40, slop_mode="middle", tokenizer_path=None, data_type="RF-HOT")
    return preprocessed_seqs

# def demo_predict_seqs(seqs, model,  progress_bar=None, print_fn=print):
#     print_fn("PREDICTING" + str(seqs.shape) + str(model) )
#     if(seqs is None):
#         print_fn("NO INPUT SEQUENCE GIVEN.")
#         return None
#     if(seqs is None or model is None):
#         print_fn("NO ML MODEL GIVEN.")
#         return None
#     y_probs = model.predict_proba(seqs)
#     y_pred  = y_probs[:, 1] if y_probs.shape[1] == 2 else y_probs[:, 0]
#     print_fn(y_probs)
#     print_fn(y_pred)
#     return y_pred

# if __name__ == "__main__":
#     app = QApplication([])
#     widget = Promotech_UI(
#         ui_path=dir_path+"/form.ui",
#         init_function=None,
#         preprocess_seqs_fn=demo_preprocess,
#         predict_seqs_fn=demo_predict_seqs,
#         predict_gen_fn=None
#     )
#     widget.show()
#     sys.exit(app.exec_())