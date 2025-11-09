"""enhance_run_model_with_pipeline_fields

Revision ID: a1b2c3d4e5f6
Revises: 83e6f8cb0bc7
Create Date: 2025-01-09 12:00:00.000000

"""

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

# revision identifiers, used by Alembic.
revision = "a1b2c3d4e5f6"
down_revision = "83e6f8cb0bc7"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # Add new columns to runs table
    op.add_column(
        "runs", sa.Column("pipeline_type", sa.String(length=50), nullable=True)
    )
    op.add_column(
        "runs", sa.Column("pipeline_template_id", sa.String(length=100), nullable=True)
    )
    op.add_column("runs", sa.Column("custom_pipeline_id", sa.Integer(), nullable=True))
    op.add_column(
        "runs",
        sa.Column(
            "parameters",
            postgresql.JSON(astext_type=sa.Text()),
            nullable=True,
            server_default="{}",
        ),
    )
    op.add_column(
        "runs",
        sa.Column(
            "input_files",
            postgresql.JSON(astext_type=sa.Text()),
            nullable=True,
            server_default="[]",
        ),
    )
    op.add_column(
        "runs",
        sa.Column(
            "output_files",
            postgresql.JSON(astext_type=sa.Text()),
            nullable=True,
            server_default="[]",
        ),
    )
    op.add_column(
        "runs",
        sa.Column(
            "input_mapping",
            postgresql.JSON(astext_type=sa.Text()),
            nullable=True,
            server_default="{}",
        ),
    )
    op.add_column("runs", sa.Column("logs", sa.Text(), nullable=True))
    op.add_column("runs", sa.Column("error_message", sa.Text(), nullable=True))
    op.add_column(
        "runs", sa.Column("progress", sa.Float(), nullable=True, server_default="0.0")
    )


def downgrade() -> None:
    # Remove added columns
    op.drop_column("runs", "progress")
    op.drop_column("runs", "error_message")
    op.drop_column("runs", "logs")
    op.drop_column("runs", "input_mapping")
    op.drop_column("runs", "output_files")
    op.drop_column("runs", "input_files")
    op.drop_column("runs", "parameters")
    op.drop_column("runs", "custom_pipeline_id")
    op.drop_column("runs", "pipeline_template_id")
    op.drop_column("runs", "pipeline_type")
