#include "myqgraphicsview.h"

MyQGraphicsView::MyQGraphicsView(QWidget *parent) :
    QGraphicsView(parent)
{
//    scene = new QGraphicsScene();
//    this->setSceneRect(50, 50, 350, 350);
//    this->setScene(scene);
//    connect(this,SIGNAL(viewClicked()),this,SLOT(saveImage()));
}

void MyQGraphicsView::mousePressEvent(QMouseEvent * e)
{
//    LOG_MSG("view mouse click start");
}

void MyQGraphicsView::mouseReleaseEvent(QMouseEvent * e)
{
//    emit viewClicked();
    if (e->button() == Qt::RightButton) {
       saveImage();
    }
}

void MyQGraphicsView::saveImage()
{
    QGraphicsScene *ascene = this->scene();
    ascene->clearSelection();                                                  // Selections would also render to the file
    ascene->setSceneRect(ascene->itemsBoundingRect());                          // Re-shrink the scene to it's bounding contents
    QImage image(ascene->sceneRect().size().toSize(), QImage::Format_ARGB32);  // Create the image with the exact size of the shrunk scene
    image.fill(Qt::transparent);                                              // Start all pixels transparent
    QPainter painter(&image);
    ascene->render(&painter);
//    char filename[] = "view.png";
//    ifield++;
//    char numstr[5];
//    sprintf(numstr,"%04d",hour);
//    for (int i=0; i<4; i++)
//        filename[11+i] = numstr[i];
    QString fileName = QFileDialog::getSaveFileName(this, tr("Image File Name"), ".", tr("Image Files (*.png)"));
    if (fileName.compare("") != 0) {
        image.save(fileName);
    }
}
