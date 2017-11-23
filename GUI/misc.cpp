#include <string>
#include <fstream>
#ifdef _WIN32
#include <windows.h>
#endif
#include <QTcpServer>
#include <QTcpSocket>
#include <QtGui>
#include <QTcpServer>
#include <QMessageBox>

#include "misc.h"
#include "log.h"
#include "transfer.h"

#include "libscell.h"

#include "global.h"

extern "C" void get_scene(int *, CELL_DATA *, bool *, double *);

LOG_USE();
char msg[2048];

class SleeperThread : public QThread
{
public:
    static void msleep(unsigned long msecs)
    {
        QThread::msleep(msecs);
    }
};

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
SocketHandler::SocketHandler(int newport, QObject *parent)
	: QThread(parent)
{
    exiting = false;
    port = newport;
	sprintf(msg,"SocketHandler: port: %d",port);
	LOG_MSG(msg);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
SocketHandler::~SocketHandler() // make sure the worker object is destroyed
{
    exiting = true;
    wait();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void SocketHandler::stop()
{
	LOG_MSG("SocketHandler::stop: set stopped");
	stopped = true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void SocketHandler::run()
{
//	QObject::moveToThread(this);
	sprintf(msg,"run: port: %d", port);
	LOG_MSG(msg);
	quint16 qport = port;
	QString addressStr = "127.0.0.1";
	QHostAddress hostAddress;
	hostAddress.setAddress(addressStr);
    tcpServer = new QTcpServer(this);
	stopped = false;
	connect(tcpServer, SIGNAL(newConnection()), this, SLOT(processor()), Qt::DirectConnection);
    if (!tcpServer->listen(hostAddress,qport)) {
 //       QMessageBox::critical(this, tr("Fortune Server"),
 //                              tr("Unable to start the server: %1.")
 //                              .arg(tcpServer->errorString()));
		sprintf(msg,"Unable to start the server: port: %d", port);
		LOG_MSG(msg);
        return;
    }
	sprintf(msg,"Listening on port: %d",tcpServer->serverPort());
	LOG_MSG(msg);
	LOG_MSG("serverAddress:");
	LOG_QMSG((tcpServer->serverAddress()).toString());
	bool timedOut = false;
	tcpServer->waitForNewConnection(-1,&timedOut);
	exec();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void SocketHandler::processor()
{
//	LOG_MSG("In processor");
    socket = tcpServer->nextPendingConnection();
	sprintf(msg,"got server connection: %p",socket);	
	LOG_MSG(msg);
    emit sh_connected();
	QString qdata;
	QByteArray ba;
	ba.resize(1024);
	while (true) {
		if (stopped) {
			LOG_MSG("SocketHandler::processor: stopped!");
			break;
		}
		socket->waitForReadyRead(100);
		int nb = socket->bytesAvailable();
		if (nb > 0) {
			ba = socket->readLine(1024);
			qdata = QString(ba);
			QStringList s = qdata.split("^",QString::SkipEmptyParts);
			for (int k=0; k<s.length(); k++) {
				emit sh_output(s[k]); // Emit signal to update GUI
				if (port == CPORT0) {
					LOG_QMSG(s[k]);
				}
			}
			if (quitMessage(qdata)) {
				sprintf(msg,"Closing connection: port: %d", port);
				LOG_MSG(msg);
		        break;
			} else {
//				LOG_MSG("No bytes yet");
			}
		}
	}
	socket->close();
	tcpServer->close();
	if (port == CPORT0) {
		emit sh_disconnected();		// Is it right that both threads emit this?
	}
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
ExecThread::ExecThread(QString infile)
{
	inputFile = infile;
}

//-----------------------------------------------------------------------------------------
// Note: the Fortran array must be SAVEd
//-----------------------------------------------------------------------------------------
void test_get_array()
{
    int i;
    int narr1, narr2;
    double *farr;
    test_array(&narr1, &narr2, &farr);
    for (i=0;i<narr1*narr2; i++)
        printf("%d %f\n",i,farr[i]);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::run()
{
	LOG_MSG("Invoking DLL...");
	int res=0;
    int hour;
	const char *infile, *outfile;
    char version[12];
	QString infile_path, outfile_path;
	int len_infile, len_outfile;
    int len_version;
    QString dll_version;
    bool cused[32];

	infile_path = inputFile;
	QString casename = QFileInfo(inputFile).baseName();
	len_infile = infile_path.length();
	std::string std_infile = infile_path.toStdString();
	infile = std_infile.c_str();
	outfile_path = casename.append(".res");
	len_outfile = outfile_path.length();
	std::string std_outfile = outfile_path.toStdString();
    outfile = std_outfile.c_str();

    get_dll_build_version(version,&len_version);
    dll_version = version;
    // We can compare the version of the actual DLL that is linked
    // with the GUI's idea of the DLL version (i.e. the version of the library the GUI was built with)
    if (dll_version != Global::DLL_build_version) {
        LOG_QMSG("Bad DLL version: " + dll_version);
        emit(badDLL(dll_version));
        return;
    }

//    test_get_array();
//    return;

	paused = false;
    LOG_MSG("call execute");
    execute(&ncpu,const_cast<char *>(infile),&len_infile,const_cast<char *>(outfile),&len_outfile,Global::blobcentre);
    LOG_MSG("did execute");
    get_dimensions(&Global::NX, &Global::NY, &Global::NZ, &nsteps, &Global::DELTA_T, &Global::MAX_CHEMO, &Global::N_EXTRA,
                   cused, &Global::dfraction, &Global::DELTA_X);
    summary_interval = int(3600./Global::DELTA_T);
    // convert DELTA_X from cm to um
    Global::DELTA_X = 10000*Global::DELTA_X;
    sprintf(msg,"exthread: nsteps: %d summary_interval: %d NX: %d DELTA_X: %f nt_VTK: %d N_EXTRA: %d",
            nsteps,summary_interval,Global::NX,Global::DELTA_X,Global::nt_vtk,Global::N_EXTRA);
    LOG_MSG(msg);
    Global::conc_nc_ex = 0;
    Global::conc_nc_ic = 0;
    hour = 0;

    mutex1.lock();
    emit setupC();
    mutex1.unlock();

    get_summary(Global::summaryData, &Global::i_hypoxia_cutoff, &Global::i_growth_cutoff);
    emit summary(hour);		// Emit signal to initialise summary plots

    for (int i=1; i <= nsteps+1; i++) {
		bool updated = false;
		if (paused && !updated) {
			snapshot();
            updated = true;
		}
        while(paused || Global::leftb) {
            sleep(100);
		}
        if (stopped) {
            res = -1;
            break;
        }

        mutex1.lock();
        simulate_step(&res);
        mutex1.unlock();
        if (res != 0) {
            LOG_MSG("simulate_step: res != 0");
            exit(1);
            break;
        }
//        LOG_MSG("did simulate_step");
        if (i%summary_interval == 0) {
            mutex1.lock();
            get_summary(Global::summaryData, &Global::i_hypoxia_cutoff, &Global::i_growth_cutoff);
            get_volprob(&Global::vol_nv, &Global::vol_v0, &Global::vol_dv, Global::volProb);
            get_oxyprob(&Global::oxy_nv, &Global::oxy_v0, &Global::oxy_dv, Global::oxyProb);
            get_concdata(&Global::conc_nvars, &Global::conc_nc_ex, &Global::conc_dx_ex, Global::concData);
            get_ic_concdata(&Global::conc_nvars, &Global::conc_nc_ic, &Global::conc_dx_ic, Global::IC_concData);
            if (Global::showingFACS || Global::recordingFACS) {
                getFACS();
            }
            mutex1.unlock();
            if (Global::showingFACS || Global::recordingFACS) {
                emit facs_update();
                emit histo_update();
            }
            hour++;
            mutex1.lock();
            emit summary(hour);		// Emit signal to update summary plots, at hourly intervals
//            summary_done.wait(&mutex3);
        }

        if (stopped) {
            res = -1;
            break;
        }
        if (i%Global::nt_vtk == 0) {
            if (Global::showingVTK || Global::recordingVTK) {
				snapshot();
                Global::istep = i;
                sleep(10);
			}
		}
        if (stopped) {
            res = -1;
            break;
        }
    }
    LOG_MSG("ExecThread::run: stopped or completed");
    stopped = true;
    snapshot();
    sleep(100);
	LOG_MSG("ExecThread::run: call terminate_run");
	terminate_run(&res);

	return;
}

/*
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::wait_to_go()
{
    for (;;) {
        sprintf(msg,"waiting");
        if (goflag || stopped) break;
    }
}
*/

//-----------------------------------------------------------------------------------------
// Note that storage for BC_list, DC_list, bond_list is provided in the GUI code
// (see mainwindow.cpp and transfer.h).  The DLL fills in the data, and the number of elements
// in the lists is returned in nBC_list, nDC_list, nbond_list.
//-----------------------------------------------------------------------------------------
void ExecThread::snapshot()
{
//    LOG_MSG("got snapshot");
    mutex2.lock();
    get_scene(&Global::ncell_list,Global::cell_list,&Global::dropped,Global::droppedcentre);
    if (Global::ncell_list > MAX_CELLS) {
        LOG_MSG("Error: MAX_CELLS exceeded");
        exit(1);
    }
    mutex2.unlock();
    emit display(); // Emit signal to update VTK display
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::saveGradient2D(int i)
{
    LOG_QMSG("saveGradient2D");
    paused = true;
    SimpleView2D *sv2D = new SimpleView2D();
    sv2D->makeFrame(i);
    paused = false;
    delete sv2D;
}

//-----------------------------------------------------------------------------------------
// Get FACS data and histogram data
//-----------------------------------------------------------------------------------------
void ExecThread::getFACS()
{
    get_nfacs(&Global::nFACS_cells);
    if (!Global::FACS_data || Global::nFACS_cells*Global::nvars_used > Global::nFACS_dim) {
        if (Global::FACS_data) free(Global::FACS_data);
        Global::nFACS_dim = 3*Global::nFACS_cells*Global::nvars_used;   // 3* to avoid excessive malloc/free
        Global::FACS_data = (double *)malloc(Global::nFACS_dim*sizeof(double));
    }
    get_facs(Global::FACS_data, Global::FACS_vmin, Global::FACS_vmax, Global::volume_scaling);

    if (!Global::histo_data || Global::nhisto_bins*Global::nvars_used > Global::nhisto_dim) {
        if (Global::histo_data) free(Global::histo_data);
        if (Global::histo_data_log) free(Global::histo_data_log);
        Global::nhisto_dim = 6*Global::nhisto_bins*Global::nvars_used;   // 2*3 to avoid excessive malloc/free (only 3* used)
        Global::histo_data = (double *)malloc(Global::nhisto_dim*sizeof(double));
        Global::histo_data_log = (double *)malloc(Global::nhisto_dim*sizeof(double));
    }

    get_histo(Global::nhisto_bins, Global::histo_data, Global::histo_vmin, Global::histo_vmax,
              Global::histo_data_log, Global::histo_vmin_log, Global::histo_vmax_log, Global::volume_scaling);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::stop()
{
	stopped = true;
	LOG_MSG("ExecThread::stop: stopped");
}

	//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::pause()
{
	paused = true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::unpause()
{
	paused = false;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool quitMessage(QString msg)
{
	if (msg.contains("__EXIT__",Qt::CaseSensitive))
		return true;
	else
		return false;
}
