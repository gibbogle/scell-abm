#ifndef SimpleView3DUI_H
#define SimpleView3DUI_H
 
#include "vtkSmartPointer.h"
#include <QMainWindow>
#include "vtkStructuredGrid.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkActor.h"
#include <QInputDialog>
#include <QMessageBox>
#include "log.h"

#ifdef __cplusplus
extern "C" {
#endif
    void get_gradient_info(int *, int *);
    void get_gradients(int *, int *, float *, int *);
#ifdef __cplusplus
}
#endif

// Forward Qt class declarations
class Ui_SimpleView3D;
 
class SimpleView3D : public QMainWindow
{
  Q_OBJECT
public:
 
  // Constructor/Destructor
  SimpleView3D();
  ~SimpleView3D() {};

  vtkSmartPointer<vtkRenderWindow> GetRenderWindow();
//  void setScale(void);
  void ShowSize(int *);
  void displayFields(void);
  void aimCamera(void);
  void create();
  void chooseParameters();
  void setParameters();
  void makeFrame(int i);

  int max_chemo;
  char msg[1024];
  int axis;
  double scale;
  int use_strength;
  int chemo_select[4];  // Currently hard-coded for 4 chemokines
  bool chemo_displayed[4];
  bool chemo_used[4];
  vtkSmartPointer<vtkRenderer> renderer;
  vtkSmartPointer<vtkRenderWindow> renWin;
  vtkRenderWindowInteractor * iren;
  vtkSmartPointer<vtkStructuredGrid> sgrid_array[4];
  vtkSmartPointer<vtkActor> sgridActor_array[4];

public slots:
 
  virtual void slotExit();
  void stateChanged_CheckBox_S1P();
  void stateChanged_CheckBox_CCL21();
  void stateChanged_CheckBox_Oxy();
  void stateChanged_CheckBox_CXCL13();
  void saveImage();

protected:
 
protected slots:
 
private:

  void CreateGradientData(vtkSmartPointer<vtkStructuredGrid> sgrid_array[], int chemo_select[], float gmax[]);
  void CreateTestData(vtkStructuredGrid* sgrid);
  void resizeEvent(QResizeEvent *);

  // Designer form
  Ui_SimpleView3D *ui;
};
 
#endif // SimpleView3DUI_H
