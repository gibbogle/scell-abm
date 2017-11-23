// myvtk.h
#ifndef MYVTK_H
#define MYVTK_H

#include <QtGui>
#include <QtCore>
#include <QIODevice>
#include <QVTKWidget.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include "vtkSphereSource.h"
#include "vtkCylinderSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
//#include <vtkMPEG2Writer.h>
#include <vtkPNGWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkBMPWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkSmartPointer.h>
#include <vtkImageCast.h>

#include <vtkAppendPolyData.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolygon.h>
#include <vtkRegularPolygonSource.h>

#include <QInputDialog>
#include <QFileDialog>
#ifdef _WIN32
#include "windows.h"
#define sleep(n) Sleep(n)   // milliseconds
#else
#include "unistd.h"
#define sleep(n) usleep(1000*n)
#endif

using namespace std;

struct cell_pos {
    int tag;
    double x, y, z;     // was int
	double diameter;
//    int state;
    int celltype;
//    int highlight;
    int status;
};
typedef cell_pos CELL_POS;

struct bond_pos {
	int BCtag;
	int DCtag;
};
typedef bond_pos BOND_POS;

struct actor_str {
    bool active;
    vtkActor *actor;
};
typedef actor_str ACTOR_TYPE;

struct colour_str {
    double r, g, b;
};
typedef colour_str COLOUR_TYPE;

#define USE_CELLTYPE_COLOUR true

class MyVTK
{
private slots:
//    void on_checkBox_CELLDISPLAY_1_toggled(bool display);
//    void on_checkBox_CELLDISPLAY_2_toggled(bool);

public:
    MyVTK(QWidget *, QWidget *);
	~MyVTK();

    void key_canvas(QWidget *);
    void createMappers();
    void get_cell_positions();
	void init();
	void cleanup();
	void unpack(int x, double *, double *, double *);
    void renderCells();
    void process_Tcells();
    void MakeWellBottom();
//    void process_Dcells();
//    void process_bonds();
	bool startPlayer(QString, QTimer *, bool);
//	bool nextFrame();
	void pause();
	void playon();
    void saveSnapshot(QString, QString, QString);
    void startRecorder(QString basename, int nframes);
    void stopRecorder();
    void recorder();
    void stop();
    void set_celltype_colour(COLOUR_TYPE *, QString str);
    void setOpacity(int celltype, int position);

    QList<CELL_POS > TCpos_list;
    QList<ACTOR_TYPE> T_Actor_list;
    vtkSmartPointer<vtkActor> sqactor;
    vtkSmartPointer<vtkPolyDataMapper> sqmapper;

    QWidget *page_VTK;
	QVTKWidget* qvtkWidget;
	vtkRenderWindow *renWin;	
	vtkRenderer* ren;
	vtkRenderWindowInteractor * iren;
    vtkPolyDataMapper *TcellMapper;
//	vtkMPEG2Writer *mpg;
//	vtkSmartPointer<vtkPNGWriter> writer;
//	vtkSmartPointer<vtkBMPWriter> writer;
    vtkSmartPointer<vtkJPEGWriter> jpgwriter;
    vtkSmartPointer<vtkPNGWriter> pngwriter;
//	vtkSmartPointer<vtkTIFFWriter> writer;
//	vtkSmartPointer<vtkImageCast> castFilter;
//	vtkWindowToImageFilter *w2img;
//	vtkSmartPointer<vtkPNGWriter> pngwriter;
//	vtkSmartPointer<vtkJPEGWriter> jpgwriter;
    vtkSmartPointer<vtkWindowToImageFilter> w2i;

	char msg[2048];
	double zoomlevel;
	double Pi;
	bool DCmotion;
	bool DCfade;
	bool first_VTK;
	bool playing;
	bool paused;
	bool save_image;
    bool display_celltype[10];
    bool dropped;
    double opacity[10];
    QColor celltype_colour[10];
    QString casename;
	int framenum;
	QTimer *timer;
	QString infile;
	QFile *playerData;
	QTextStream *playerStream;

    // Image recorder
//    bool record;
//    QString record_basename;
//    int record_nframes;
//    int record_it;
//    QTemporaryFile * tempFile;

public slots:
//    void on_checkBox_CELLDISPLAY_1_toggled(bool display);
//    void on_checkBox_CELLDISPLAY_2_toggled(bool);

};

#endif
