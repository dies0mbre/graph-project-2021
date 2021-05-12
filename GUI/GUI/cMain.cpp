
//************************************************************************
//************************************************************************
// BertImage
// created: february 2006
// updated: may 2011
// How to load, display, process and save images with wxWindows
// Author: Pascal Bertolino, GIPSA-lab laboratory, Grenoble - France
// Email pascal.bertolino@gipsa-lab.fr
// Web http://www.gipsa-lab.inpg.fr/~pascal.bertolino/
// tested with xWidget 2.8.7 under Linux and Windows
//************************************************************************
//************************************************************************

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include <wx/image.h>
#include <wx/file.h>
#include <wx/bitmap.h>

#define APP_NAME "Image Segmentation"

enum
{
	ID_QUIT = 1,
	ID_ABOUT,
	ID_LOAD,
	ID_SAVE,
	ID_PROCESS,
	ID_BEST_SIZE
};

//************************************************************************
//************************************************************************
// Canvas class (where we display the image)
//************************************************************************
//************************************************************************

//------------------------------------------------------------------------
class MyCanvas : public wxPanel
	//------------------------------------------------------------------------
{
public:
	MyCanvas(wxWindow* parent, wxWindowID, const wxPoint& pos, const wxSize& size);
	~MyCanvas();
	void LoadImage(wxString fileName);
	void SaveImage(wxString fileName);
	void ProcessImage();
	void BestSize();

private:
	int m_imageWidth;
	int m_imageHeight;
	wxBitmap m_imageBitmap;	// used to display the image
	wxImage* m_imageRGB;		// used to load the image
	unsigned char* m_myImage;	// used to process the image

	void OnPaint(wxPaintEvent& event);

	DECLARE_EVENT_TABLE()
};

BEGIN_EVENT_TABLE(MyCanvas, wxPanel)
EVT_PAINT(MyCanvas::OnPaint)
END_EVENT_TABLE()

//------------------------------------------------------------------------
MyCanvas::MyCanvas(wxWindow* parent, wxWindowID id,
	const wxPoint& pos, const wxSize& size)
	: wxPanel(parent, id, pos, size, wxSUNKEN_BORDER)
	//------------------------------------------------------------------------
{
	m_myImage = NULL;
	m_imageRGB = NULL;
}

//------------------------------------------------------------------------
MyCanvas::~MyCanvas()
//------------------------------------------------------------------------
{
	if (m_myImage)
		free(m_myImage);
	if (m_imageRGB)
		delete m_imageRGB;
}

//------------------------------------------------------------------------
void MyCanvas::LoadImage(wxString fileName)
//------------------------------------------------------------------------
{
	if (m_myImage)
		free(m_myImage);
	if (m_imageRGB)
		delete m_imageRGB;

	// open image dialog box
	m_imageRGB = new wxImage(fileName, wxBITMAP_TYPE_ANY, -1); // ANY => can load many image formats
	m_imageBitmap = wxBitmap(*m_imageRGB, -1); // ...to get the corresponding bitmap

	m_imageWidth = m_imageRGB->GetWidth();
	m_imageHeight = m_imageRGB->GetHeight();

	m_myImage = (unsigned char*)malloc(m_imageWidth * m_imageHeight * 3);
	memcpy(m_myImage, m_imageRGB->GetData(), m_imageWidth * m_imageHeight * 3);

	// update GUI size
	SetSize(m_imageWidth, m_imageHeight);
	GetParent()->SetClientSize(GetSize());

	// update display
	Refresh(false);
}

//------------------------------------------------------------------------
void MyCanvas::SaveImage(wxString fileName)
//------------------------------------------------------------------------
{
	bool b;

	wxImage* tempImage = new wxImage(m_imageWidth, m_imageHeight, m_myImage, true); // lend my image buffer...
	b = tempImage->SaveFile(fileName);
	delete(tempImage);		// buffer not needed any more

	if (!b)
		wxMessageBox(wxT("A problem occured during saving"));
}

//------------------------------------------------------------------------
void MyCanvas::ProcessImage()
//------------------------------------------------------------------------
// example of fast and trivial process (negative)
// you can replace it with your own
// you can also use methods from the wxImage class itself
{
	long int i = m_imageWidth * m_imageHeight * 3;

	// m_myImage is a monodimentional vector of pixels (RGBRGB...)
	while (i--)
		m_myImage[i] = 255 - m_myImage[i];

	Refresh(false); // update display
}

//------------------------------------------------------------------------
void MyCanvas::BestSize()
//------------------------------------------------------------------------
{
	SetSize(m_imageWidth, m_imageHeight);	// ideal size for canvas
	GetParent()->SetClientSize(GetSize());	// force the main frame to show the whole canvas
}

//------------------------------------------------------------------------
void MyCanvas::OnPaint(wxPaintEvent& WXUNUSED(event))
//------------------------------------------------------------------------
// update the main window content
{
	wxImage* tempImage;  // the bridge between my image buffer and the bitmap to display

	wxPaintDC dc(this);

	if (m_myImage)
	{
		tempImage = new wxImage(m_imageWidth, m_imageHeight, m_myImage, true); // lend my image buffer...
		m_imageBitmap = wxBitmap(*tempImage, -1); // ...to get the corresponding bitmap
		delete(tempImage);		// buffer not needed any more
		dc.DrawBitmap(m_imageBitmap, 0, 0);
	}
}

//************************************************************************
//************************************************************************
// Frame class (the main window)
//************************************************************************
//************************************************************************

//------------------------------------------------------------------------
class MyFrame : public wxFrame
	//------------------------------------------------------------------------
{
public:
	MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size);

	// Event handlers
protected:
	void OnQuit(wxCommandEvent& event);
	void OnAbout(wxCommandEvent& event);
	void OnOpenImage(wxCommandEvent& WXUNUSED(event));
	void OnSaveImage(wxCommandEvent& WXUNUSED(event));
	void OnProcessImage(wxCommandEvent& WXUNUSED(event));
	void OnClose(wxCloseEvent& event);
	void OnBestSize(wxCommandEvent& WXUNUSED(event));

	MyCanvas* m_canvas; // the canvas inside the main frame
	bool m_imageLoaded;
	DECLARE_EVENT_TABLE()
};


BEGIN_EVENT_TABLE(MyFrame, wxFrame)
EVT_MENU(ID_LOAD, MyFrame::OnOpenImage)
EVT_MENU(ID_SAVE, MyFrame::OnSaveImage)
EVT_MENU(ID_PROCESS, MyFrame::OnProcessImage)
EVT_MENU(ID_BEST_SIZE, MyFrame::OnBestSize)
EVT_MENU(ID_QUIT, MyFrame::OnQuit)
EVT_MENU(ID_ABOUT, MyFrame::OnAbout)
EVT_CLOSE(MyFrame::OnClose)

END_EVENT_TABLE()

//------------------------------------------------------------------------
MyFrame::MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size)
	: wxFrame((wxFrame*)NULL, -1, title, pos, size)
	//------------------------------------------------------------------------
{
	wxMenu* file_menu = new wxMenu();
	file_menu->Append(ID_LOAD, _T("&Open image..."));
	file_menu->Append(ID_PROCESS, _T("&Process image"));
	file_menu->Append(ID_SAVE, _T("&Save image as..."));
	file_menu->Append(ID_BEST_SIZE, _T("&Best size"));
	file_menu->AppendSeparator();
	file_menu->Append(ID_ABOUT, _T("&About..."));
	file_menu->AppendSeparator();
	file_menu->Append(ID_QUIT, _T("&Exit"));

	wxMenuBar* menuBar = new wxMenuBar();
	menuBar->Append(file_menu, _T("&File"));
	SetMenuBar(menuBar);

	// create the canvas that will manage the image
	m_canvas = new MyCanvas(this, -1, wxDefaultPosition, wxDefaultSize);
	m_imageLoaded = false;
	Centre();
}


//------------------------------------------------------------------------
void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	Close(true);
}

//------------------------------------------------------------------------
void MyFrame::OnClose(wxCloseEvent& event)
//------------------------------------------------------------------------
{
	delete m_canvas;
	event.Skip();
}

//------------------------------------------------------------------------
void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	wxMessageBox(_T("How to \n\n- load\n- display\n- process\n- save\n\nan image with wxWidgets (2.8.7)\n\nPascal Bertolino - GIPSA-lab, Grenoble - France\npascal.bertolino@gipsa-lab.fr"),
		_T(APP_NAME), wxOK | wxICON_INFORMATION);
}

//------------------------------------------------------------------------
void MyFrame::OnProcessImage(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	if (m_imageLoaded)
		m_canvas->ProcessImage();
}

//------------------------------------------------------------------------
void MyFrame::OnOpenImage(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	wxBitmap bitmap;

	wxString filename = wxFileSelector(_T("Select file"), _T(""), _T(""), _T(""), _T("All files (*.*)|*.*"));
	if (!filename.empty())
	{
		m_canvas->LoadImage(filename);
		m_imageLoaded = true;
	}
}

//------------------------------------------------------------------------
void MyFrame::OnSaveImage(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	//	char str[128] = "" ; // proposed file name

	if (!m_imageLoaded)
		return;

	wxString filename = wxFileSelector(_T("Save image as"), _T(""), _T(""), _T("*.bmp"), _T("BMP files (*.bmp)|*.bmp|GIF files (*gif)|*.gif|JPEG files (*jpg)|*.jpg|PNG files (*png)|*.png|TIFF files (*tif)|*.tif|XPM files (*xpm)|*.xpm|All files (*.*)|*.*"), wxFD_SAVE);
	if (!filename.empty())
		m_canvas->SaveImage(filename);
}

//------------------------------------------------------------------------
void MyFrame::OnBestSize(wxCommandEvent& WXUNUSED(event))
//------------------------------------------------------------------------
{
	m_canvas->BestSize();
}

//************************************************************************
//************************************************************************
// Application class
//************************************************************************
//************************************************************************

//------------------------------------------------------------------------
class MyApp : public wxApp
	//------------------------------------------------------------------------
{
	virtual bool OnInit();
};

IMPLEMENT_APP(MyApp) // macro that contains the main() function


//------------------------------------------------------------------------
bool MyApp::OnInit()
//------------------------------------------------------------------------
{
	//support all available image formats
	wxInitAllImageHandlers();

	MyFrame* frame = new MyFrame(_T(APP_NAME), wxDefaultPosition, wxSize(400, 300));
	frame->Show(true);
	SetTopWindow(frame);
	return true;
}
