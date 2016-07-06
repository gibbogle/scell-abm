#include "qdebug.h"
#include "field.h"
#include "log.h"

#include "global.h"

LOG_USE();

#ifdef __DISPLAY768
#define CANVAS_WIDTH 640
#else
#define CANVAS_WIDTH 696
#endif

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
Field::Field(QWidget *aParent) : QWidget(aParent)
{
    field_page = aParent;
    axis = Z_AXIS;
    fraction = 0;
    cell_constituent = OXYGEN;
    field_constituent = OXYGEN;
    slice_changed = true;
    setConcPlot(false);
    setVolPlot(false);
    setOxyPlot(false);
    pGconc = NULL;
    pGvol = NULL;
    pGoxy = NULL;
    show_cells = true;
    ifield = 0;
    view = new MyQGraphicsView(field_page);
    scene = new QGraphicsScene(QRect(0, 0, CANVAS_WIDTH, CANVAS_WIDTH));
    vbox_cell_constituent = NULL;
    vbox_field_constituent = NULL;
    vbox_cell_max_concentration = NULL;
    buttonGroup_cell_constituent = new QButtonGroup;
    buttonGroup_field_constituent = new QButtonGroup;
    line_maxConc_list.clear();
    cell_constituent_rb_list.clear();
    field_constituent_rb_list.clear();
//    data = NULL;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
Field::~Field()
{
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool Field::isConcPlot()
{
    return useConcPlot;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::setConcPlot(bool status)
{
    useConcPlot = status;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool Field::isVolPlot()
{
    return useVolPlot;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::setVolPlot(bool status)
{
    useVolPlot = status;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool Field::isOxyPlot()
{
    return useOxyPlot;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::setOxyPlot(bool status)
{
    useOxyPlot = status;
}

//------------------------------------------------------------------------------------------------
// To create the group of radiobuttons for cell constituent selection.
// This uses information about active constituents fetched from the DLL.
// Need to distinguish field constituent from cell constituent, because for the cell we are
// interested in CFSE, volume, O2byvol in addition to the dissolved constituents:
// oxygen, glucose, drugA, drugB, metabolites...
//------------------------------------------------------------------------------------------------
void Field::setCellConstituentButtons(QGroupBox *gbox, QButtonGroup *bg, QVBoxLayout **vbox, QList<QRadioButton *> *rb_list, QString tag)
{
    int ivar;
    QString name, str, ivar_str;
    QRadioButton *rb;

    LOG_QMSG("setCellConstituentButtons: " + tag);
    if (rb_list->length() != 0) {
        LOG_MSG("rb_list not NULL, delete it");
        for (ivar=0; ivar<rb_list->length(); ivar++) {
            rb = (*rb_list)[ivar];
            bg->removeButton(rb);
            delete rb;
        }
        rb_list->clear();
    }
    if (!*vbox) {
        LOG_MSG("vbox = NULL, create it");
        *vbox = new QVBoxLayout;
        gbox->setLayout(*vbox);
    }
    name = "rb_cell_constituent_"+tag;
    LOG_QMSG(name);
    sprintf(msg,"rb_list: %p vbox: %p bg: %p nvars_used: %d",rb_list,*vbox,bg,Global::nvars_used);
    LOG_MSG(msg);
    for (ivar=0; ivar<Global::nvars_used; ivar++) {
        ivar_str = QString::number(ivar);
        str = Global::var_string[ivar];
        rb = new QRadioButton;
        rb->setText(str);
        rb->setObjectName(name+ivar_str);
        (*vbox)->addWidget(rb);
        rb->setEnabled(true);
        bg->addButton(rb,ivar);
        rb_list->append(rb);
//        LOG_QMSG(rb->objectName());
    }
    LOG_MSG("added buttons");
    if (tag.contains("FACS")) {
        (*rb_list)[0]->setChecked(true);   // CFSE
    } else {
        (*rb_list)[1]->setChecked(true);   // Oxygen
    }
    QRect rect = gbox->geometry();
    rect.setHeight(25*Global::nvars_used);
    gbox->setGeometry(rect);
    gbox->show();
}


//------------------------------------------------------------------------------------------------
// To create the group of radiobuttons for field constituent selection.
// This uses information about active constituents fetched from the DLL.
// Need to distinguish field constituent from cell constituent, because for the
// field we have only the dissolved constituents:
// oxygen, glucose, drugA, drugB, metabolites...
//------------------------------------------------------------------------------------------------
void Field::setFieldConstituentButtons(QGroupBox *gbox, QButtonGroup *bg, QVBoxLayout **vbox, QList<QRadioButton *> *rb_list, QString tag)
{
    int ivar;
    QString name, str, ivar_str;
    QRadioButton *rb;

    Global::nfieldvars_used = Global::nvars_used - Global::N_EXTRA;
//    LOG_QMSG("setFieldConstituentButtons: " + tag + " nfieldvars_used: "
//             + QString::number(Global::nfieldvars_used));
    if (rb_list->length() != 0) {
        LOG_MSG("rb_list not NULL, delete it");
        for (ivar=0; ivar<rb_list->length(); ivar++) {
            rb = (*rb_list)[ivar];
            bg->removeButton(rb);
            delete rb;
        }
        rb_list->clear();
    }
    if (!*vbox) {
        LOG_MSG("vbox = NULL, create it");
        *vbox = new QVBoxLayout;
        gbox->setLayout(*vbox);
    }
    name = "rb_field_constituent_"+tag;
//    LOG_QMSG(name);
    for (ivar=0; ivar<Global::nfieldvars_used; ivar++) {
        ivar_str = QString::number(ivar);
        str = Global::var_string[ivar+1];
        rb = new QRadioButton;
        rb->setText(str);
        QString objname = name+ivar_str;
        rb->setObjectName(objname);
        (*vbox)->addWidget(rb);
        rb->setEnabled(true);
        bg->addButton(rb,ivar);
        rb_list->append(rb);
//        QRadioButton *chkrb = gbox->findChild<QRadioButton *>(objname);
//        if (chkrb) {
//            QString chkstr = chkrb->objectName();
//            LOG_QMSG(chkstr);
//        } else {
//            chkrb = (*rb_list)[ivar];
//            LOG_QMSG("findChild failed, but: " + chkrb->objectName());
//        }
    }
    (*rb_list)[0]->setChecked(true);   // Oxygen
    QRect rect = gbox->geometry();
    rect.setHeight(25*(Global::nfieldvars_used + 1));
    gbox->setGeometry(rect);
    gbox->show();
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setMaxConcentrations(QGroupBox *gbox)
{
    int ivar;
    QLineEdit *line_maxConc;

    if (!vbox_cell_max_concentration) {
        LOG_MSG("vbox_cell_max_concentration = NULL, create it");
        vbox_cell_max_concentration = new QVBoxLayout;
        gbox->setLayout(vbox_cell_max_concentration);
    }

    if (line_maxConc_list.length() != 0) {
        for (ivar=0; ivar < line_maxConc_list.length(); ivar++) {
            line_maxConc = line_maxConc_list[ivar];
            vbox_cell_max_concentration->removeWidget(line_maxConc);
            line_maxConc->setVisible(false);
        }
    }
    line_maxConc_list.clear();

    for (ivar=0; ivar<Global::nvars_used; ivar++) {
        line_maxConc = new QLineEdit;
        line_maxConc->setObjectName("maxConc"+ivar);
        line_maxConc->setText("1.0");
        line_maxConc->setMaximumWidth(50);
        line_maxConc->setEnabled(true);
        vbox_cell_max_concentration->addWidget(line_maxConc);
        line_maxConc_list.append(line_maxConc);
    }

    line_maxConc_list[OXYGEN]->setText("0.18");
    line_maxConc_list[GLUCOSE]->setText("9.0");
    QRect rect = gbox->geometry();
    rect.setHeight(25*Global::nvars_used);
    gbox->setGeometry(rect);
    gbox->show();
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::selectCellConstituent()
{
    int iconst, res;
    QStringList items;

//    LOG_MSG("selectCellConstituent");
    for (iconst=0; iconst<Global::nvars_used; iconst++) {
        if (iconst == cell_constituent) continue;
        items << Global::var_string[iconst];
    }
    bool ok;
    QString item = QInputDialog::getItem(this, tr("QInputDialog::getItem()"),
                                          tr("Constituent:"), items, 0, false, &ok);
    if (ok && !item.isEmpty()) {
        for (iconst=0; iconst<Global::nvars_used; iconst++) {
            if (item == Global::var_string[iconst]) {
                cell_constituent = iconst;
//                LOG_QMSG("selectCellConstituent: " + QString::number(cell_constituent));
            }
        }
    }
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setCellConstituent(QAbstractButton *button)
{
    int res;
    LOG_MSG("setCellConstituent");
    int prev_constituent = cell_constituent;
    QString text = button->text();
    for (int ivar=0; ivar<Global::nvars_used; ivar++) {
        if (text == Global::var_string[ivar]) {
            cell_constituent = ivar;
        }
    }

    if (cell_constituent != prev_constituent) {
        LOG_MSG("setCellConstituent -> displayField");
        slice_changed = true;
        displayField(hour,&res);
        LOG_MSG("did displayField");
    }
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::selectFieldConstituent()
{
    int iconst, res;
    QStringList items;

//    LOG_MSG("selectFieldConstituent");
    for (iconst=1; iconst<Global::nfieldvars_used+1; iconst++) {
        if (iconst == field_constituent) continue;
        items << Global::var_string[iconst];
    }
    bool ok;
    QString item = QInputDialog::getItem(this, tr("QInputDialog::getItem()"),
                                          tr("Constituent:"), items, 0, false, &ok);
    if (ok && !item.isEmpty()) {
        for (iconst=1; iconst<Global::nfieldvars_used+1; iconst++) {
            if (item == Global::var_string[iconst]) {
                field_constituent = iconst;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setFieldConstituent(QAbstractButton *button)
{
    int res;
//    LOG_MSG("setFieldConstituent");
    int prev_constituent = field_constituent;
    QString text = button->text();
    for (int ivar=0; ivar<Global::nfieldvars_used; ivar++) {
        if (text == Global::var_string[ivar+1]) {
            field_constituent = ivar+1;
        }
    }

    if (field_constituent != prev_constituent) {
//        LOG_MSG("setFieldConstituent");
        slice_changed = true;
        displayField(hour,&res);
    }
}


//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setPlane(QAbstractButton *button)
{
    int res;
    LOG_MSG("setPlane");
    QString text = button->text();
	int prev_axis = axis;
    if (text.compare("X-Y") == 0)
        axis = Z_AXIS;
    else if (text.compare("Y-Z") == 0)
        axis = X_AXIS;
    else if (text.compare("X-Z") == 0)
        axis = Y_AXIS;
	if (axis != prev_axis) {
        slice_changed = true;
        LOG_MSG("setPlane");
        displayField(hour,&res);
	}
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setFraction(QString text)
{
    int res;
	double prev_fraction = fraction;
	fraction = text.toDouble();
	if (fraction != prev_fraction) {
        displayField(hour,&res);
	}
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setSliceChanged()
{
    slice_changed = true;
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setSaveImages(bool save)
{
    save_images = save;
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setUseLogScale(bool use_logscale)
{
    use_log = use_logscale;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::chooseFieldColor(double c, double cmin, double cmax, bool use_logscale, int rgbcol[])
{
    double f, denom, logcmin, logc;
    int rgb_lo[3], rgb_hi[3], i;
    bool copy_scell = false;

    if (copy_scell) {
        rgbcol[0] = 0;
        if (cmax > 0) {
            rgbcol[1] = 255*min(c,cmax)/cmax;
            rgbcol[2] = 255*min(c,cmax)/cmax;
        } else {
            rgbcol[1] = 0;
            rgbcol[2] = 0;
        }
        return;
    }
    if (use_logscale) {
        if (cmin == cmax) {
            f = 1;
        } else {
//            if (cmin > 0.0001)
//                logcmin = log(cmin);
//            else
//                logcmin = 1.0e6;
//            if (c > 0.0001)
//                logc = log(c);
//            else
//                logc = 1.0e6;
            cmin = max(cmin, 0.00001);
            c = max(c, 0.00001);
            denom = (log(cmax) - log(cmin));
            if (denom < 0.001)
                f = 1;
            else
                f = (log(c) - log(cmin))/denom;
        }
    } else {
        f = c/cmax;
    }
//    if (cell_constituent == OXYGEN) {
    if (field_constituent == OXYGEN) {
        rgb_hi[0] =   0; rgb_hi[1] =   0; rgb_hi[2] = 0;
        rgb_lo[0] = 255; rgb_lo[1] =   0; rgb_lo[2] = 0;
        for (i=0; i<3; i++) {
            rgbcol[i] = int((1-f)*rgb_lo[i] + f*rgb_hi[i]);
            if (rgbcol[i] < 0 || rgbcol[i] > 255) {
                sprintf(msg,"chooseFieldColor: %f %f %f %f %d %d",c,cmin,cmax,f,i,rgbcol[i]);
                LOG_MSG(msg);
                exit(1);
            }
        }
    } else {
        rgb_hi[0] =   0; rgb_hi[1] =   255; rgb_hi[2] = 255;
        rgb_lo[0] =   0; rgb_lo[1] =   0; rgb_lo[2] = 0;
        for (i=0; i<3; i++) {
            rgbcol[i] = int((1-f)*rgb_lo[i] + f*rgb_hi[i]);
            if (rgbcol[i] < 0 || rgbcol[i] > 255) {
                sprintf(msg,"chooseFieldColor: %f %f %f %f %d %d",c,cmin,cmax,f,i,rgbcol[i]);
                LOG_MSG(msg);
                exit(1);
            }
        }
    }
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::chooseRateColor(double f, int rgbcol[])
{
    int rgb_lo[3], rgb_hi[3], i;

    rgb_hi[0] = 0; rgb_hi[1] = 255; rgb_hi[2] = 0;
    rgb_lo[0] = 0; rgb_lo[1] =  64; rgb_lo[2] = 0;
    for (i=0; i<3; i++) {
        rgbcol[i] = int((1-f)*rgb_lo[i] + f*rgb_hi[i]);
    }
}

//-----------------------------------------------------------------------------------------
// Now 'constituent' is an index of the active constituents: 0 - nvars_used
//-----------------------------------------------------------------------------------------
void Field::getTitle(int iconst, QString *title)
{
    QString name = Global::var_string[iconst];
    if (Global::GUI_to_DLL_index[iconst] <= Global::MAX_CHEMO) {
        *title = name + " Concentration";
    } else {
        *title = name;
    }
}

//-----------------------------------------------------------------------------------------
// New version, site/cell size is fixed, the blob grows
//-----------------------------------------------------------------------------------------
void Field::displayField(int hr, int *res)
{
//    QGraphicsScene* scene = new QGraphicsScene(QRect(0, 0, CANVAS_WIDTH, CANVAS_WIDTH));
    QBrush brush;
    int i, k, ix, iy, iz, w, rgbcol[3], ichemo, ixyz;
    double xp, yp, x, y, d, C, cmin, cmax, rmax, valmin;
    double a, b, Wc, dx, Wx, radius;
    int Nc, NX, NY, NZ;
    double beta = 1.0;

    ichemo = Global::GUI_to_DLL_index[field_constituent];
//    LOG_QMSG("displayField: field: " + QString::number(field_constituent) + " --> " + QString::number(ichemo));
    use_log = false;    // temporary
    *res = 0;
    if (hr >= 0) hour = hr;
	if (slice_changed) {
        get_fielddata(&axis, &fraction, &fdata, &ixyz, res);
        if (*res != 0) {
            LOG_MSG("Error: get_fielddata: FAILED");
            return;
        }
//        sprintf(msg,"fdata: %d %d %d %d ncells: %d",fdata.NX,fdata.NY,fdata.NZ,fdata.NCONST,fdata.ncells);
//        LOG_MSG(msg);
        slice_changed = false;
    }

    scene->clear();

    cmin = 1.0e10;
    cmax = 0;
    rmax = 1;

    NX = fdata.NX;
    NY = fdata.NY;
    NZ = fdata.NZ;
    dx = fdata.DX;
    /*
    if (axis == X_AXIS) {           // Y-Z plane
        nx = NY;
        ny = NZ;
        ix = (NX+1)/2 - 1;
        for (iy=0; iy<NY; iy++) {
            for (iz=0; iz<NZ; iz++) {
                k = (ichemo-1)*NZ*NY*NX + iz*NY*NX + iy*NX + ix;    // index(ix,iy,iz,ichemo);
                C = fdata.Cave[k];
//                cmin = MIN(MAX(cmin,0),C);
                cmax = MAX(cmax,C);
                cmin = 0;
//                cmax = line_maxConc_list[cell_constituent]->text().toDouble();
            }
        }
    } else if (axis == Y_AXIS) {   // X-Z plane
        nx = NX;
        ny = NZ;
        iy = (NY+1)/2 - 1;
    } else if (axis == Z_AXIS) {   // X-Y plane
        nx = NX;
        ny = NY;
        iz = (NZ+1)/2 - 1;
    }
    */

    Wx = (NX-1)*dx;
    Wc = CANVAS_WIDTH;
    a = Wc/(beta*Wx);
    b = Wc/2 - a*Wx/2;
//    sprintf(msg,"Wx: %f Wc: %f a: %f b: %f",Wx,Wc,a,b);
//    LOG_MSG(msg);
    brush.setStyle(Qt::SolidPattern);
    brush.setColor(QColor(0,0,0));
    scene->addRect(0,0,CANVAS_WIDTH,CANVAS_WIDTH,Qt::NoPen, brush);
    view->setScene(scene);
    view->setGeometry(QRect(0, 0, CANVAS_WIDTH+4, CANVAS_WIDTH+4));

//    sprintf(msg,"displayField: field constituent: %d ichemo: %d cmax: %f",field_constituent,ichemo,cmax);
//    LOG_MSG(msg);
    cmax = line_maxConc_list[field_constituent]->text().toDouble();
//    sprintf(msg,"displayField: field constituent: %d ichemo: %d cmax: %f",field_constituent,ichemo,cmax);
//    LOG_MSG(msg);
    rgbcol[0] = 0;
    rgbcol[1] = 0;
    rgbcol[2] = 0;
    w = a*dx + 0.5;
    valmin = 1.0e10;
    if (axis == X_AXIS) {           // Y-Z plane
        ix = ixyz;
        for (iy=0; iy<NY; iy++) {
            x = iy*dx;
            for (iz=0; iz<NZ; iz++) {
                y = (NZ-1-iz)*dx;
                k = (ichemo-1)*NZ*NY*NX + iz*NY*NX + iy*NX + ix;    // index(ix,iy,iz,ichemo);
                C = fdata.Cave[k];
                valmin = MIN(C,valmin);
                rgbcol[1] = 255*min(C,cmax)/cmax;
                rgbcol[2] = 255*min(C,cmax)/cmax;
                xp = int(a*x + b);
                yp = int(a*y + b);
    //        chooseFieldColor(data[i].conc[ichemo],cmin,cmax,use_log,rgbcol);
    //        sprintf(msg,"c: %f %f %f rgbcol: %d %d %d",data[i].conc[constituent],cmin,cmax,rgbcol[0],rgbcol[1],rgbcol[2]);
    //        LOG_MSG(msg);
                brush.setColor(QColor(rgbcol[0],rgbcol[1],rgbcol[2]));
                scene->addRect(xp,yp,w,w,Qt::NoPen, brush);
            }
        }
    } else if (axis == Y_AXIS) {           // X-Z plane
        iy = ixyz;
        for (ix=0; ix<NX; ix++) {
            x = ix*dx;
            for (iz=0; iz<NZ; iz++) {
                y = (NZ-1-iz)*dx;
                k = (ichemo-1)*NZ*NY*NX + iz*NY*NX + iy*NX + ix;    // index(ix,iy,iz,ichemo);
                C = fdata.Cave[k];
                valmin = MIN(C,valmin);
                rgbcol[1] = 255*min(C,cmax)/cmax;
                rgbcol[2] = 255*min(C,cmax)/cmax;
                xp = int(a*x + b);
                yp = int(a*y + b);
                brush.setColor(QColor(rgbcol[0],rgbcol[1],rgbcol[2]));
                scene->addRect(xp,yp,w,w,Qt::NoPen, brush);
            }
        }
    } else if (axis == Z_AXIS) {           // X-Y plane
        iz = ixyz;
        for (ix=0; ix<NX; ix++) {
            x = ix*dx;
            for (iy=0; iy<NY; iy++) {
                y = iy*dx;
                k = (ichemo-1)*NZ*NY*NX + iz*NY*NX + iy*NX + ix;    // index(ix,iy,iz,ichemo);
                C = fdata.Cave[k];
//                if (ix == NX/2) {
//                    sprintf(msg,"%d %d %d %d %8.3f",ichemo,ix,iy,iz,C);
//                    LOG_MSG(msg);
//                }
                valmin = MIN(C,valmin);
                rgbcol[1] = 255*min(C,cmax)/cmax;
                rgbcol[2] = 255*min(C,cmax)/cmax;
                xp = int(a*x + b);
                yp = int(a*y + b);
                brush.setColor(QColor(rgbcol[0],rgbcol[1],rgbcol[2]));
                scene->addRect(xp,yp,w,w,Qt::NoPen, brush);
            }
        }
    }
//    sprintf(msg,"axis: %d valmin: %f",axis,valmin);
//    LOG_MSG(msg);

    if (!show_cells) return;
//    ichemo = Global::GUI_to_DLL_index[cell_constituent];
//    LOG_QMSG("displayField: cell: " + QString::number(cell_constituent) + " --> " + QString::number(ichemo));
//    LOG_QMSG("displayField: nc: " + QString::number(fdata.ncells));
    for (i=0; i<fdata.ncells; i++) {
        x = fdata.cell_data[i].centre[0];
        y = fdata.cell_data[i].centre[1];
        radius = fdata.cell_data[i].radius;
        xp = a*x + b;
        yp = a*y + b;
        d = 2*a*radius;
//        sprintf(msg,"Cell: %d x,y: %f %f radius: %f xp,yp: %f %f",i,x,y,radius,xp,yp);
//        LOG_MSG(msg);
        if (fdata.cell_data[i].status == 0) {
            rgbcol[0] = 0;
            rgbcol[1] = 200;
            rgbcol[2] = 32;
        } else if (fdata.cell_data[i].status == 1) {
            rgbcol[0] = 50;
            rgbcol[1] = 100;
            rgbcol[2] = 32;
        } else if (fdata.cell_data[i].status == 2 || fdata.cell_data[i].status == 4) {
            rgbcol[0] = 0;
            rgbcol[1] = 0;
            rgbcol[2] = 255;
        } else if (fdata.cell_data[i].status == 3) {
            rgbcol[0] = 255;
            rgbcol[1] = 0;
            rgbcol[2] = 255;
        } else if (fdata.cell_data[i].status >= 10) {   // tagged to die of treatment
            rgbcol[0] = 255;
            rgbcol[1] = 150;
            rgbcol[2] = 0;
        }
        brush.setColor(QColor(rgbcol[0],rgbcol[1],rgbcol[2]));
        scene->addEllipse(xp,yp,d,d,Qt::NoPen, brush);
    }

    double w_scalebar = a*0.01 + b; // 100um = 0.01cm
    double scalebar0 = a*1*dx + b;
    QPen pen;
    QFont font;
    pen.setBrush(Qt::black);
    pen.setWidth(3);
    scene->addLine(scalebar0, scalebar0, scalebar0+w_scalebar, scalebar0, pen);
    QGraphicsTextItem *scalebar_text = scene->addText("100 um",font);
    scalebar_text->setPos(scalebar0,1.4*scalebar0);
    view->show();
    return;

    if (save_images) {
        scene->clearSelection();                                                  // Selections would also render to the file
        scene->setSceneRect(scene->itemsBoundingRect());                          // Re-shrink the scene to it's bounding contents
        QImage image(scene->sceneRect().size().toSize(), QImage::Format_ARGB32);  // Create the image with the exact size of the shrunk scene
        image.fill(Qt::transparent);                                              // Start all pixels transparent

        QPainter painter(&image);
        scene->render(&painter);
        ifield++;
        char filename[] = "image/field0000.png";
        char numstr[5];
        sprintf(numstr,"%04d",hour);
        for (int i=0; i<4; i++)
            filename[11+i] = numstr[i];
        image.save(filename);
    }
}


/*
//-----------------------------------------------------------------------------------------
// New version, site/cell size is fixed, the blob grows
//-----------------------------------------------------------------------------------------
void Field::displayField(int hr, int *res)
{
    QGraphicsScene* scene = new QGraphicsScene(QRect(0, 0, CANVAS_WIDTH, CANVAS_WIDTH));
    QBrush brush;
    int i, xindex, yindex, ix, iy, w, rgbcol[3], ichemo;
    double xp, yp, d0, d, volume, scale, cmin, cmax, rmax;
    double a, b, Wc;
    int Nc;

    ichemo = Global::GUI_to_DLL_index[field_constituent];
//    LOG_QMSG("displayField: " + QString::number(field_constituent) + "-->" + QString::number(ichemo));
    use_log = false;    // temporary
    *res = 0;
    hour = hr;
	if (slice_changed) {
        get_fieldinfo(&NX, &axis, &fraction, &nsites, &nconst, const_used, res);    // Note: const_used[] is no longer used
        if (*res != 0) {
            printf("Error: get_fieldinfo: FAILED\n");
            LOG_MSG("Error: get_fieldinfo: FAILED");
            return;
        }
        if (nconst != MAX_CONC) {
            sprintf(msg,"Error: displayField: get_fieldinfo: MAX_CONC != MAX_CHEMO: %d %d",MAX_CONC,Global::MAX_CHEMO);
            LOG_MSG(msg);
            exit(1);
        }
        if (NEXTRA != Global::N_EXTRA) {
            sprintf(msg,"Error: displayField: NEXTRA != N_EXTRA: %d %d",NEXTRA,Global::N_EXTRA);
            LOG_MSG(msg);
            exit(1);
        }
        if (this->data) {
            free(this->data);
        }
        this->data = (FIELD_DATA *)malloc(nsites*sizeof(FIELD_DATA));
        get_fielddata(&axis, &fraction, &nsites, &nconst, this->data, res);
        if (*res != 0) {
            LOG_MSG("Error: get_fielddata: FAILED");
            return;
        }
        slice_changed = false;
    }
    if (axis == X_AXIS) {           // Y-Z plane
        xindex = 1;
        yindex = 2;
    } else if (axis == Y_AXIS) {   // X-Z plane
        xindex = 0;
        yindex = 2;
    } else if (axis == Z_AXIS) {   // X-Y plane
        xindex = 0;
        yindex = 1;
    }


// NX = size of lattice
// Nc = # of sites to fill the canvas from side to side (or top to bottom) = (2/3)NX
// Wc = canvas width (pixels)
// w = site width = Wc/Nc
// xp = a.ix + b
// yp = a.iy + b
// blob centre at (NX/2,NX/2) maps to canvas centre at (Wc/2,Wc/2)
// => Wc/2 = a.NX/2 + b
// The width of Nc sites maps to the canvas width
// => Wc = a.Nc
// => a = Wc/Nc, b = Wc/2 - a.NX/2
    Nc = (2*NX)/3;
    Wc = CANVAS_WIDTH;
    w = Wc/Nc;
    a = w;
    b = Wc/2 - a*NX/2;
    d0 = w*Global::dfraction;
    cmin = 1.0e10;
    cmax = 0;
    rmax = 1;
    for (i=0; i<nsites; i++) {
        cmin = MIN(MAX(cmin,0),data[i].conc[ichemo]);
        cmax = MAX(cmax,data[i].conc[ichemo]);
        // Flip it
//        data[i].site[yindex] = NX - data[i].site[yindex];
    }
    brush.setStyle(Qt::SolidPattern);
    brush.setColor(QColor(0,0,0));
    scene->addRect(0,0,CANVAS_WIDTH,CANVAS_WIDTH,Qt::NoPen, brush);
    view->setScene(scene);
    view->setGeometry(QRect(0, 0, 700, 700));
    if (cmax == 0) {
        sprintf(msg,"displayField: ichemo: %d cmax=0",ichemo);
        LOG_MSG(msg);
    } else {
        for (i=0; i<nsites; i++) {
            ix = this->data[i].site[xindex];
            iy = this->data[i].site[yindex];
            xp = int(a*ix + b - w);
            yp = int(a*iy + b - w);
            chooseFieldColor(data[i].conc[ichemo],cmin,cmax,use_log,rgbcol);
    //        sprintf(msg,"c: %f %f %f rgbcol: %d %d %d",data[i].conc[cell_constituent],cmin,cmax,rgbcol[0],rgbcol[1],rgbcol[2]);
    //        LOG_MSG(msg);
    //        rgbcol[1] = 0;
            brush.setColor(QColor(rgbcol[0],rgbcol[1],rgbcol[2]));
            scene->addRect(xp,yp,w,w,Qt::NoPen, brush);
        }
    }
    for (i=0; i<nsites; i++) {
        ix = this->data[i].site[xindex];
        iy = this->data[i].site[yindex];
        xp = int(a*ix + b - w);
        yp = int(a*iy + b - w);
        volume = this->data[i].volume;      // = 0 if there is no cell
        if (volume > 0) {
            double f = 0;
            scale = pow(volume,0.3333);
            d = scale*d0;   // fix this - need to change d0
            if (this->data[i].state == 0) {
                f = data[i].conc[GROWTH_RATE]/rmax;      // currently colour cells only by growth rate, now fraction of max growth rate
                chooseRateColor(f,rgbcol);
            } else if (this->data[i].state == 1) {      // hypoxic (O2 < hypoxic threshold)
                rgbcol[0] = 50;
                rgbcol[1] = 100;
                rgbcol[2] = 32;
            } else if (this->data[i].state == 2) {      // anoxia_tag
                rgbcol[0] = 0;
                rgbcol[1] = 0;
                rgbcol[2] = 255;
            } else if (this->data[i].state == 3) {      // mitosis
                rgbcol[0] = 255;
                rgbcol[1] = 0;
                rgbcol[2] = 255;
            } else if (this->data[i].state >= 10) {     // tagged to die from radiation or drug
                rgbcol[0] = 255;
                rgbcol[1] = 150;
                rgbcol[2] = 0;
            }
            brush.setColor(QColor(rgbcol[0],rgbcol[1],rgbcol[2]));
            scene->addEllipse(xp+(w-d)/2,yp+(w-d)/2,d,d,Qt::NoPen, brush);
//            sprintf(msg,"i: %6d f: %6.3f vol: %6.3f",i,f,volume);
//            LOG_MSG(msg);
        }
    }
    view->show();
    if (save_images) {
        scene->clearSelection();                                                  // Selections would also render to the file
        scene->setSceneRect(scene->itemsBoundingRect());                          // Re-shrink the scene to it's bounding contents
        QImage image(scene->sceneRect().size().toSize(), QImage::Format_ARGB32);  // Create the image with the exact size of the shrunk scene
        image.fill(Qt::transparent);                                              // Start all pixels transparent

        QPainter painter(&image);
        scene->render(&painter);
        ifield++;
        char filename[] = "image/field0000.png";
        char numstr[5];
        sprintf(numstr,"%04d",hour);
        for (int i=0; i<4; i++)
            filename[11+i] = numstr[i];
        LOG_QMSG("image.save");
        LOG_MSG(filename);
        image.save(filename);
    }
}
*/

