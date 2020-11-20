/**
 * line intersection test program
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date 11-Nov-2020
 * @license: see 'LICENSE' file
 */

#ifndef __LINES_GUI_H__
#define __LINES_GUI_H__


#include <QMainWindow>
#include <QLabel>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsItem>

#include <memory>
#include <vector>

#include "geo_algos.h"
using t_real = double;


class Vertex : public QGraphicsItem
{
public:
	Vertex(const QPointF& pos, double rad = 15.);
	virtual ~Vertex();

	virtual QRectF boundingRect() const override;
	virtual void paint(QPainter*, const QStyleOptionGraphicsItem*, QWidget*) override;

private:
	double m_rad = 15.;
};


enum class IntersectionCalculationMethod
{
	DIRECT,
	SWEEP,
};


class LinesView : public QGraphicsView
{Q_OBJECT
public:
	using t_vec = m::vec<t_real, std::vector>;
	using t_mat = m::mat<t_real, std::vector>;

	static const constexpr t_real g_eps = 1e-5;


public:
	LinesView(QGraphicsScene *scene=nullptr, QWidget *parent=nullptr);
	virtual ~LinesView();

	LinesView(LinesView&) = delete;
	const LinesView& operator=(const LinesView&) const = delete;

	void AddVertex(const QPointF& pos);
	void ClearVertices();
	const std::vector<Vertex*>& GetVertexElems() const { return m_elems_vertices; }

	void UpdateAll();
	void UpdateLines();
	void UpdateIntersections();

	void SetIntersectionCalculationMethod(IntersectionCalculationMethod m);

	void SaveVoro(const std::string& filename) const;

protected:
	virtual void mousePressEvent(QMouseEvent *evt) override;
	virtual void mouseReleaseEvent(QMouseEvent *evt) override;
	virtual void mouseMoveEvent(QMouseEvent *evt) override;

	virtual void resizeEvent(QResizeEvent *evt) override;

private:
	QGraphicsScene *m_scene = nullptr;

	std::vector<Vertex*> m_elems_vertices{};
	std::vector<QGraphicsItem*> m_elems_lines{}, m_elems_inters{};

	bool m_dragging = false;

	std::vector<std::pair<t_vec, t_vec>> m_lines{};

	IntersectionCalculationMethod m_intersectioncalculationmethod = IntersectionCalculationMethod::SWEEP;

signals:
	void SignalMouseCoordinates(double x, double y);
};



class LinesWnd : public QMainWindow
{
public:
	using QMainWindow::QMainWindow;

	LinesWnd(QWidget* pParent = nullptr);
	~LinesWnd();

	void SetStatusMessage(const QString& msg);

private:
	virtual void closeEvent(QCloseEvent *) override;

private:
	std::shared_ptr<QGraphicsScene> m_scene;
	std::shared_ptr<LinesView> m_view;
	std::shared_ptr<QLabel> m_statusLabel;
};


#endif
